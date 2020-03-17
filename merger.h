//
// Created by cothrax on 3/14/20.
//

#ifndef LIBJPEG_MERGER_H
#define LIBJPEG_MERGER_H

#include "param.h"
#include "refseq.h"
#include<vector>
using std::vector;
const int MAX_READ_DUP = 3000;

struct Hit_Can
{

    Hit_Can()
    {
        chr = 0;
        _end = 0;
        _begin = 0;
    }

    // Index of ref
    ref_id_t chr;

    // Location of alignment end
    ref_loc_t _end;

    // Location of alignment begin
    ref_loc_t _begin;

};

class Merger
{
private:
    struct Bucket {ref_loc_t l, r; bool exist, valid; };

    int _seed_size;
    int _read_length;
    int _ref_num;
    bit64_t node_num;
    int recruit_length;

    vector<ref_loc_t> *hits;
//    ref_loc_t *l, *r, *cnt;
    vector<ref_loc_t> stack;
    Bucket *buckets;
    vector<Hit_Can> *_perfect_hit;
    RefSeq *_ref;

public:
    Merger();
    ~Merger();
    void init(int seed_size, int ref_num, bit64_t max_ref_length);
    inline void put_into_bucket(ref_loc_t loc, ref_loc_t hit, bool if_jud);
    inline void checkout(ref_loc_t loc, vector<ref_loc_t> &v);
    inline void expand(ref_loc_t loc, ref_id_t ref_id);
    void collect(ref_id_t ref_id);
    void run(vector<Hit_Can> &perfect_hit, RefSeq &ref);

    void push_back(ref_id_t ref_id, ref_loc_t hit)
    {
//        printf("push_back: %d %d\n", ref_id, hit);
        hits[ref_id].push_back(hit);
    }
    void set_read_length(int read_length)
    {
        _read_length = read_length;
        recruit_length = _read_length + (_seed_size << 1);
    }
};


#endif //LIBJPEG_MERGER_H
