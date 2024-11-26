#include <nanobind/nanobind.h>
#include <iostream>
#include "simdjson.h"
#include <string>
#include <nanobind/stl/string.h>
#include <vector>
#include <parallel_hashmap/phmap.h>
#include <nanobind/stl/vector.h>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <fstream>
#include <zlib.h>
#include <kseq++/seqio.hpp>
#include <hashutil.hpp>
#include <atomic>
#include <thread>
#include <deque>
#include <condition_variable>
#include <mutex>

using namespace klibpp;
using json = nlohmann::json;

using phmap::parallel_flat_hash_map;

namespace nb = nanobind;
using namespace std;
using namespace simdjson;
using namespace nb::literals;

bool valid_file(std::string path)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (fp == NULL)
    {
        return false;
    }
    fclose(fp);
    return true;
}

class Dehasher
{

private:
    int kSize;
    vector<string> sig_paths;
    phmap::parallel_flat_hash_map<
        uint64_t, string,
        std::hash<uint64_t>,
        std::equal_to<uint64_t>,
        std::allocator<std::pair<const uint64_t, string>>,
        /* Number of submaps */ 6,
        /* Mutex type */ std::mutex>
        hash_to_kmer;

    phmap::flat_hash_set<uint64_t> all_hashes;

    // This to save all hashes in all sigs
    void fetch_hashes()
    {
        for (const string &sig_path : sig_paths)
        {
            vector<uint64_t> hashes;
            get_hashes_from_sig(sig_path, hashes);
            for (uint64_t hash_val : hashes)
            {
                this->all_hashes.insert(hash_val); // Now works with flat_hash_set
            }
        }
        cout << "Total hashes: " << this->all_hashes.size() << endl;
    }

    void get_hashes_from_large_sig(std::string sigpath, vector<uint64_t> &hashes)
    {

        std::ifstream sig_stream(sigpath);
        json data = json::parse(sig_stream);

        for (auto const &sketch : data)
        {
            for (auto const &signature : sketch["signatures"])
            {
                int _ksize = (int)signature["ksize"];
                if (this->kSize == _ksize)
                {
                    for (auto const &hash_val : signature["mins"])
                    {
                        hashes.emplace_back(hash_val);
                    }
                }
            }
        }
    }

    void get_hashes_from_sig(std::string sig_path, vector<uint64_t> &hashes)
    {
        ondemand::parser parser;
        padded_string json = padded_string::load(sig_path);
        ondemand::document sig = parser.iterate(json);

        ondemand::array sketches = sig.get_array();

        for (auto sketch : sketches)
        {
            ondemand::array signatures = sketch.find_field("signatures");

            for (auto signature : signatures)
            {
                int _ksize = (int)signature["ksize"].get_int64();
                if (this->kSize == _ksize)
                {
                    ondemand::array mins = signature["mins"].get_array();
                    for (uint64_t hash_val : mins)
                    {
                        hashes.emplace_back(hash_val);
                    }
                }
            }
        }
    }

public:
    Dehasher(int kSize, vector<string> sig_paths)
    {
        this->kSize = kSize;
        for (string sig_path : sig_paths)
        {
            if (valid_file(sig_path))
            {
                this->sig_paths.push_back(sig_path);
            }
            else
            {
                throw invalid_argument("Invalid file path: " + sig_path);
            }
        }
        this->fetch_hashes();
    }

    void map_kmer_to_hashes_single_fasta_parallel(std::string fasta_path, int num_threads)
    {
        auto *murmurHasher = new MumurHasher(42);
        std::atomic<bool> all_hashes_found(false);
        std::atomic<bool> done_reading(false);
        std::atomic<size_t> hashes_found(0);
        std::atomic<size_t> total_reads_processed(0);    // Atomic counter for total reads processed
        std::atomic<size_t> next_progress_update(10000); // Next milestone for progress update

        std::mutex progress_mutex; // Mutex to synchronize progress updates

        // Initialize queues and synchronization primitives
        std::vector<std::deque<std::string>> thread_queues(num_threads);
        std::vector<std::mutex> queue_mutexes(num_threads);
        std::vector<std::condition_variable> queue_cvs(num_threads);

        // Launch worker threads
        std::vector<std::thread> threads;
        for (int t = 0; t < num_threads; ++t)
        {
            threads.emplace_back([&, t]()
                                 {
                size_t thread_reads_processed = 0;
                while (!all_hashes_found)
                {
                    std::string seq;
                    {
                        std::unique_lock<std::mutex> lock(queue_mutexes[t]);
                        queue_cvs[t].wait(lock, [&]() {
                            return !thread_queues[t].empty() || all_hashes_found || done_reading;
                        });
                        if (all_hashes_found || (thread_queues[t].empty() && done_reading))
                            break;
                        seq = std::move(thread_queues[t].front());
                        thread_queues[t].pop_front();
                    }

                    const size_t seq_length = seq.size();
                    const char* seq_data = seq.c_str();

                    for (size_t i = 0; i <= seq_length - this->kSize && !all_hashes_found; ++i)
                    {
                        const char* kmer_ptr = seq_data + i;
                        std::string kmer(kmer_ptr, this->kSize);
                        uint64_t hash = murmurHasher->hash(kmer);

                        if (this->all_hashes.contains(hash))
                        {
                            // Insert into hash_to_kmer using try_emplace
                            auto [iterator, inserted] = this->hash_to_kmer.try_emplace(hash, kmer);

                            if (inserted)
                            {
                                size_t current_hashes_found = ++hashes_found;
                                if (current_hashes_found == this->all_hashes.size())
                                {
                                    all_hashes_found = true;
                                    break;
                                }
                            }
                        }
                    }

                    ++thread_reads_processed;
                    size_t total_processed = ++total_reads_processed;

                    // Update progress every 10,000 reads
                    if (total_processed >= next_progress_update)
                    {
                        std::lock_guard<std::mutex> progress_lock(progress_mutex);
                        if (total_processed >= next_progress_update)
                        {
                            std::cout << "\rReads processed: " << total_processed << std::flush;
                            next_progress_update += 10000;
                        }
                    }
                } });
        }

        std::cout << "Processing " << fasta_path << "..." << std::endl;

        // Main thread reads sequences and assigns them to threads
        KSeq record;
        SeqStreamIn iss(fasta_path.c_str());
        size_t seq_index = 0;
        while (iss >> record && !all_hashes_found)
        {
            int thread_idx = seq_index % num_threads;
            {
                std::unique_lock<std::mutex> lock(queue_mutexes[thread_idx]);
                thread_queues[thread_idx].emplace_back(std::move(record.seq));
            }
            queue_cvs[thread_idx].notify_one();
            ++seq_index;
        }

        // Indicate that reading is done
        done_reading = true;

        // Notify all threads to finish
        for (int t = 0; t < num_threads; ++t)
        {
            queue_cvs[t].notify_one();
        }

        // Wait for threads to finish
        for (auto &thread : threads)
        {
            thread.join();
        }

        // Final progress update
        std::cout << "\rReads processed: " << total_reads_processed << std::endl;
    }

    void map_kmer_to_hashes_multi_fasta_parallel(vector<string> fasta_paths, int num_threads)
    {
        auto *murmurHasher = new MumurHasher(42);
        std::atomic<bool> all_hashes_found(false);
        std::atomic<size_t> hashes_found(0);
        std::atomic<size_t> total_reads_processed(0);
        std::atomic<size_t> next_progress_update(10000);
        std::mutex progress_mutex;

        // Determine the number of threads to use
        size_t num_files = fasta_paths.size();
        size_t num_threads_to_use = std::min(num_threads, static_cast<int>(num_files));

        std::atomic<size_t> next_file_index(0);

        // Launch worker threads
        std::vector<std::thread> threads;
        for (size_t t = 0; t < num_threads_to_use; ++t)
        {
            threads.emplace_back([&, t]()
                                 {
                while (!all_hashes_found)
                {
                    size_t file_index = next_file_index++;
                    if (file_index >= fasta_paths.size())
                        break;

                    const std::string& fasta_path = fasta_paths[file_index];
                    KSeq record;
                    SeqStreamIn iss(fasta_path.c_str());
                    while (iss >> record && !all_hashes_found)
                    {
                        const size_t seq_length = record.seq.size();
                        const char* seq_data = record.seq.c_str();

                        for (size_t i = 0; i <= seq_length - this->kSize && !all_hashes_found; ++i)
                        {
                            const char* kmer_ptr = seq_data + i;
                            std::string kmer(kmer_ptr, this->kSize);
                            uint64_t hash = murmurHasher->hash(kmer);

                            if (this->all_hashes.contains(hash))
                            {
                                // Insert into hash_to_kmer using try_emplace
                                auto [iterator, inserted] = this->hash_to_kmer.try_emplace(hash, kmer);

                                if (inserted)
                                {
                                    size_t current_hashes_found = ++hashes_found;
                                    if (current_hashes_found == this->all_hashes.size())
                                    {
                                        all_hashes_found = true;
                                        break;
                                    }
                                }
                            }
                        }

                        size_t total_processed = ++total_reads_processed;

                        // Update progress every 10,000 reads
                        if (total_processed >= next_progress_update)
                        {
                            std::lock_guard<std::mutex> progress_lock(progress_mutex);
                            if (total_processed >= next_progress_update)
                            {
                                std::cout << "\rReads processed: " << total_processed << std::flush;
                                next_progress_update += 10000;
                            }
                        }
                    }
                } });
        }

        std::cout << "Processing multiple FASTA files..." << std::endl;

        // Wait for all threads to finish
        for (auto &thread : threads)
        {
            thread.join();
        }

        // Final progress update
        std::cout << "\rReads processed: " << total_reads_processed << std::endl;
    }

    void map_kmer_to_hashes_single_fasta(std::string fasta_path)
    {
        auto *murmurHasher = new MumurHasher(42);

        size_t total_reads_processed = 0;
        bool all_hashes_found = false;

        KSeq record;
        SeqStreamIn iss(fasta_path.c_str());
        while (iss >> record && !all_hashes_found)
        {
            for (size_t i = 0; i <= record.seq.size() - this->kSize && !all_hashes_found; i++)
            {
                const char *kmer_ptr = &record.seq[i];
                uint64_t hash = murmurHasher->hash(kmer_ptr, this->kSize);
                if (this->all_hashes.find(hash) != this->all_hashes.end())
                {
                    this->hash_to_kmer[hash] = std::string(kmer_ptr, this->kSize);
                    if (this->hash_to_kmer.size() == this->all_hashes.size())
                    {
                        all_hashes_found = true;
                    }
                }
            }
            total_reads_processed++;
        }
        std::cout << "\rReads processed: ~" << total_reads_processed << std::flush;
        std::cout << std::endl;
    }

    unordered_map<uint64_t, string> get_hash_to_kmer()
    {
        unordered_map<uint64_t, string> hash_to_kmer;
        for (auto const &hash_kmer : this->hash_to_kmer)
        {
            hash_to_kmer[hash_kmer.first] = hash_kmer.second;
        }
        return hash_to_kmer;
    }

    void dump_kmers_to_file(string file_path)
    {
        ofstream out_file(file_path);
        for (auto const &hash_kmer : this->hash_to_kmer)
        {
            out_file << hash_kmer.first << "\t" << hash_kmer.second << endl;
        }
        out_file.close();
    }
};

NB_MODULE(_dehasher_impl, m)
{
    nb::class_<Dehasher>(m, "Dehasher")
        .def(nb::init<int, vector<string>>(), "kSize"_a, "sig_paths"_a)
        .def("map_kmer_to_hashes_single_fasta_parallel", &Dehasher::map_kmer_to_hashes_single_fasta_parallel, "fasta_path"_a, "num_threads"_a)
        .def("map_kmer_to_hashes_multi_fasta_parallel", &Dehasher::map_kmer_to_hashes_multi_fasta_parallel, "fasta_paths"_a, "num_threads"_a)
        .def("map_kmer_to_hashes_single_fasta", &Dehasher::map_kmer_to_hashes_single_fasta, "fasta_path"_a)
        .def("get_hash_to_kmer", &Dehasher::get_hash_to_kmer)
        .def("dump_kmers_to_file", &Dehasher::dump_kmers_to_file, "file_path"_a);
}
