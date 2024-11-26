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
#include <unordered_map>
#include <zlib.h>
#include <kseq++/seqio.hpp>
#include <hashutil.hpp>

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
    int chunk_size;
    phmap::flat_hash_map<uint64_t, std::string> hash_to_kmer;
    phmap::flat_hash_map<uint64_t, bool> all_hashes;

    // This to save all hashes in all sigs
    void fetch_hashes()
    {
        for (string sig_path : sig_paths)
        {
            vector<uint64_t> hashes;
            get_hashes_from_sig(sig_path, hashes);
            for (uint64_t hash_val : hashes)
            {
                this->all_hashes[hash_val] = true;
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
    Dehasher(int kSize, vector<string> sig_paths, int chunk_size = 1)
    {
        this->kSize = kSize;
        this->chunk_size = chunk_size;
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

    void map_kmer_to_hashes_multi_fasta(vector<string> fasta_paths)
    {

        // initialize the hasher
        auto *murmurHasher = new MumurHasher(42);

        for (string fasta_path : fasta_paths)
        {
            // initialize kseq
            KSeq record;
            // fasta path must be C string
            SeqStreamIn iss(fasta_path.c_str());
            while (iss >> record)
            {
                for (int i = 0; i < record.seq.size() - this->kSize + 1; i++)
                {
                    std::string kmer = record.seq.substr(i, this->kSize);
                    uint64_t hash = murmurHasher->hash(kmer);
                    if (this->all_hashes.find(hash) != this->all_hashes.end())
                    {
                        this->hash_to_kmer[hash] = kmer;
                    }
                }
            }
        }
    }

    void map_kmer_to_hashes_single_fasta(std::string fasta_path)
    {
        auto *murmurHasher = new MumurHasher(42);

        size_t total_reads_processed = 0;

        KSeq record;
        SeqStreamIn iss(fasta_path.c_str());
        while (iss >> record)
        {
            for (int i = 0; i < record.seq.size() - this->kSize + 1; i++)
            {
                std::string kmer = record.seq.substr(i, this->kSize);
                uint64_t hash = murmurHasher->hash(kmer);
                if (this->all_hashes.find(hash) != this->all_hashes.end())
                {
                    this->hash_to_kmer[hash] = kmer;
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
        .def(nb::init<int, vector<string>, int>(), "kSize"_a, "sig_paths"_a, "chunk_size"_a = 1)
        .def("map_kmer_to_hashes_single_fasta", &Dehasher::map_kmer_to_hashes_single_fasta, "fasta_path"_a)
        .def("map_kmer_to_hashes_multi_fasta", &Dehasher::map_kmer_to_hashes_multi_fasta, "fasta_paths"_a)
        .def("get_hash_to_kmer", &Dehasher::get_hash_to_kmer)
        .def("dump_kmers_to_file", &Dehasher::dump_kmers_to_file, "file_path"_a);
}
