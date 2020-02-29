//
// Created by cothrax on 2/24/20.
//

#ifndef LIBJPEG_STAT_H
#define LIBJPEG_STAT_H

#include<fstream>
#include<iostream>
#include<cstdio>
#include<cstring>
#include<vector>
#include "param.h"
#include "refseq.h"

using namespace std;
typedef vector<pair<bit64_t, bit64_t> > StatVector;
const char default_hits_fn[20] = "hits.dat";
const char default_matches_fn[20] = "matches.dat";

bit64_t transform(bit64_t read_id, bit64_t read_loc, bit64_t if_reversed = 0);
void inv_transform(bit64_t val, ref_id_t &read_id, ref_loc_t &read_loc, bool &if_reversed);
void comp_reverse(char *seq);


class StatIO
{
protected:
    void write(ofstream &fout, StatVector &vec, bit64_t size);
    ofstream fout_hit, fout_match;

public:
//    vector<bit64_t> read_info;
    bit64_t _if_reversed, _read_id;
    StatVector hit_vec, match_vec;
    vector<vector<int> > hit_history;
    int *sort_index;
    vector<int> *Ahit;

    StatIO(const char *hits_fn = default_hits_fn, const char *matches_fn = default_matches_fn);
    ~StatIO();
    void flush(bit64_t batch_size = 100);
    bool comp(int x, int y);
    void argsort(vector<bit64_t> &Ahit);

    void set_read_info(bit64_t read_id, bit64_t if_reversed)
    {
        _read_id = read_id;
        _if_reversed = if_reversed;
    }
    void push_back_hit(bit64_t read_loc, bit64_t ref_id, bit64_t ref_loc)
    {
        hit_vec.push_back(make_pair(transform(_read_id, read_loc, _if_reversed), transform(ref_id, ref_loc)));
    }
    void push_back_match(int cand_id)
    {
        for(int i = 0; i < hit_history[cand_id].size(); i++) match_vec.push_back(hit_vec[hit_history[cand_id][i]]);
    }
    void push_back_vector()
    {
        hit_history.push_back(vector<int>());
    }
    void push_into_vector(int val)
     {
        hit_history[hit_history.size()-1].push_back(sort_index[val]);
    }
};


class DataHolder
{
public:
    bit64_t *len, *ofs;
    char *seq;
    int _seed_size;
    bit64_t seq_num;

    explicit DataHolder(int seed_size=11): _seed_size(seed_size) {};
    ~DataHolder();
    void get_len(const char *filename, int size);
    void load(const char *filename, int size = -1);
    void get(char *res, ref_id_t id, ref_loc_t loc, bool if_reversed);
    void get2(char *res, ref_id_t id, ref_loc_t start, ref_loc_t end);
};


class StatConverter
{
public:
    ifstream fin;
    explicit StatConverter(const char *filename);
    ~StatConverter();
    bool next_pair(bit64_t &read_val, bit64_t &ref_val);
    void convert(DataHolder &ref, DataHolder &reads, const char *output_fn);
};

#endif //LIBJPEG_STAT_H
