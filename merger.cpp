//
// Created by cothrax on 3/14/20.
//

#include "merger.h"
#include<cstring>

//void Merger::put_into_bucket(ref_loc_t loc, ref_loc_t hit)
//{
//    if(buckets[loc].cnt++ == 0)
//    {
//        buckets[loc].l = buckets[loc].r = hit;
//        stack.push_back(loc);
//    }
//    else
//    {
//        if(hit < buckets[loc].l) buckets[loc].l = hit;
//        if(hit > buckets[loc].r) buckets[loc].r = hit;
//    }
//}

inline void Merger::put_into_bucket(ref_loc_t loc, ref_loc_t hit, bool if_jud)
{
    if(!buckets[loc].exist)
    {
        buckets[loc].exist = true;
        if(if_jud) buckets[loc].valid = buckets[loc+1].exist || buckets[loc-1].exist;
        buckets[loc].l = buckets[loc].r = hit;
        stack.push_back(loc);
    }
    else
    {
        if(if_jud) buckets[loc].valid = true;
        if(hit < buckets[loc].l) buckets[loc].l = hit;
        if(hit > buckets[loc].r) buckets[loc].r = hit;
    }
}

void Merger::collect(ref_id_t ref_id)
{
    int i;
    vector<ref_loc_t> &v = hits[ref_id];
    if(!v.size()) return;

    for(i = 0; i < v.size(); i++) put_into_bucket(v[i] / _seed_size + 1, v[i], true);
    v.clear();
    for(i = 0; i < stack.size(); i++) checkout(stack[i], v);
//        if(checkout(stack[i]))
//        {
//            v.push_back(buckets[stack[i]].l);
//            v.push_back(buckets[stack[i]].r);
//        }
    stack.clear();
    for(i = 0; i < v.size(); i++) put_into_bucket(v[i] / recruit_length + 1, v[i], false);
    for(i = 0; i < stack.size(); i++) expand(stack[i], ref_id);
    stack.clear();
    v.resize(0);
}

void Merger::run(vector<Hit_Can> &perfect_hit, RefSeq &ref)
{
    _perfect_hit = &perfect_hit;
    _ref = &ref;
    for(int i = 0; i < _ref_num; i++) collect(i);
    stack.resize(0);
}

inline void Merger::checkout(ref_loc_t loc, vector<ref_loc_t> &v)
{
    if(buckets[loc].valid) v.push_back(buckets[loc].l);
    buckets[loc].valid = buckets[loc].exist = false;

//    if(buckets[loc].cnt > 1)
//    {
//        buckets[loc].cnt = 0;
//        return true;
//    }
//    if(buckets[loc].cnt == 1)
//    {
//        if(buckets[loc+1].cnt && buckets[loc+1].l < _seed_size + buckets[loc].r)
//        {
//            buckets[loc].cnt = buckets[loc+1].cnt = 0;
//            return true;
//        }
//        if(buckets[loc-1].cnt && buckets[loc].l < _seed_size + buckets[loc-1].r)
//        {
//            buckets[loc].cnt = buckets[loc-1].cnt = 0;
//            return true;
//        }
//    }
//    buckets[loc].cnt = 0;
//    return false;
}

//inline void Merger::expand(ref_loc_t loc, ref_id_t ref_id)
//{
//    if(!buckets[loc].cnt) return;
//    ref_loc_t l_loc = loc, r_loc = loc;
//
//    buckets[loc].cnt = 0;
//    while(buckets[l_loc-1].cnt && buckets[l_loc].l < recruit_length + buckets[l_loc-1].r) buckets[--l_loc].cnt = 0;
//    while(buckets[r_loc+1].cnt && buckets[r_loc+1].l < recruit_length + buckets[r_loc].r) buckets[++r_loc].cnt = 0;
//    Hit_Can can;
//    can.chr = ref_id;
//    can._begin = buckets[l_loc].l > _seed_size ? buckets[l_loc].l - _seed_size : 0;
//    can._end = buckets[r_loc].r + _read_length + (_seed_size << 1);
//    if(can._end - can._begin > MAX_READ_DUP) can._end = can._begin + MAX_READ_DUP;
//    if(can._end >= _ref->title[ref_id].size) can._end = _ref->title[ref_id].size;
//
//    _perfect_hit->push_back(can);
//}

inline void Merger::expand(ref_loc_t loc, ref_id_t ref_id)
{
    if(!buckets[loc].exist) return;
    ref_loc_t l_loc = loc, r_loc = loc;

    buckets[loc].exist = 0;
    while(buckets[l_loc-1].exist && buckets[l_loc].l < recruit_length + buckets[l_loc-1].r) buckets[--l_loc].exist = false;
    while(buckets[r_loc+1].exist && buckets[r_loc+1].l < recruit_length + buckets[r_loc].r) buckets[++r_loc].exist = false;
    Hit_Can can;
    can.chr = ref_id;
    can._begin = buckets[l_loc].l > _seed_size ? buckets[l_loc].l - _seed_size : 0;
    can._end = buckets[r_loc].r + _read_length + (_seed_size << 1);
    if(can._end - can._begin > MAX_READ_DUP) can._end = can._begin + MAX_READ_DUP;
    if(can._end >= _ref->title[ref_id].size) can._end = _ref->title[ref_id].size;

    _perfect_hit->push_back(can);
}

void Merger::init(int seed_size, int ref_num, bit64_t max_ref_length)
{
    _seed_size = seed_size;
//    _read_length = read_length;
    _ref_num = ref_num;
    node_num = max_ref_length / seed_size + 2;

    hits = new vector<ref_loc_t>[ref_num];
    for(int i = 0; i < ref_num; i++) hits[i] = vector<ref_loc_t>();
    buckets = new Bucket[node_num];
//    buckets++;
    memset(buckets, 0, sizeof(Bucket) * node_num);
    stack = vector<ref_loc_t>();
}

Merger::Merger()
{
    hits = 0;
    buckets = 0;
    _perfect_hit = 0;
}

Merger::~Merger()
{
    delete [] hits;
    delete [] buckets;
}
