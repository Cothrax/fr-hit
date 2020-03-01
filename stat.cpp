//
// Created by cothrax on 2/24/20.
//

#include "stat.h"
#include<algorithm>
#include<cmath>
using namespace std;
//char _seq[700000000];

// --------------------------------- helper func ---------------------------------
bit64_t transform(bit64_t read_id, bit64_t read_loc, bit64_t if_reversed)
{
//    if(((if_reversed << 63) | (read_id << 31) | read_loc) == 384400929203)
//        cout << "!" << read_id << " " << read_loc << " " << if_reversed << endl;
    return (if_reversed << 63) | (read_id << 32) | read_loc;
}

void inv_transform(bit64_t val, ref_id_t &read_id, ref_loc_t &read_loc, bool &if_reversed)
{
    static bit64_t one = 1;
    if_reversed = val>>63;
    val &= ~(one << 63);
    read_loc = val & ((one<<32)-1);
    read_id = (val>>32) & ((one<<32)-1);
}

void comp_reverse(char *seq)
{
    int len = strlen(seq);
    for(int i = 0; i < len; i++)
        switch(seq[i])
        {
            case 'A': seq[i] = 'T'; break;
            case 'G': seq[i] = 'C'; break;
            case 'C': seq[i] = 'G'; break;
            case 'T': seq[i] = 'A';
        }
    for(int i = 0; i < len / 2; i++) swap(seq[i], seq[len-1-i]);
}


// --------------------------------- StatIO ---------------------------------
StatIO::StatIO(const char *hits_fn, const char *matches_fn)
{
    fout_hit = ofstream(hits_fn, ios_base::out | ios_base::binary);
    fout_match = ofstream(matches_fn, ios_base::out | ios_base::binary);
    sort_index = nullptr;
}

StatIO::~StatIO()
{
    fout_hit.close();
    fout_match.close();
}

void StatIO::write(ofstream &fout, StatVector &vec, bit64_t size)
{
    for(int i = 0; i < size; i++)
    {
        fout.write((char *) &vec[i].first, sizeof(vec[i].first));
        fout.write((char *) &vec[i].second, sizeof(vec[i].second));
    }
    fout.flush();
}

void StatIO::flush(bit64_t batch_size)
{
    write(fout_hit, hit_vec, batch_size);
    write(fout_match, match_vec, match_vec.size());

//    int ofs = 0;
//    for(int j = 0; j < match_size.size(); j++)
//    {
//        fout.write((char *) &match_size[i], sizeof(match_size[i]));
//        for(int i = 0; i < match_size[j]; i++)
//        {
//            fout.write((char *) &match_vec[i].first, sizeof(match_vec[i].first));
//            fout.write((char *) &match_vec[i].second, sizeof(match_vec[i].second));
//        }
//    }
//    for(int i = 0; i < match_vec.size(); i++) cout << "match: " << match_vec[i].first << " " << match_vec[i].second << endl;

    for(int i = 0; i < hit_history.size(); i++) hit_history[i].resize(0);
    hit_history.resize(0);
    match_vec.resize(0);
    hit_vec.resize(0);
    delete [] sort_index;
}

bool StatIO::comp(int x, int y) { return (*Ahit)[x] < (*Ahit)[y]; }

void StatIO::argsort(vector <bit64_t> &Ahit)
{
    sort_index = new int[Ahit.size()];
    pair<bit64_t, int> *tmp_arr = new pair<bit64_t, int>[Ahit.size()];
    for(int i = 0; i < Ahit.size(); i++) tmp_arr[i] = make_pair(Ahit[i], i);
    sort(tmp_arr, tmp_arr+Ahit.size());
    for(int i = 0; i < Ahit.size(); i++) sort_index[i] = tmp_arr[i].second;
}

// --------------------------------- DataHolder ---------------------------------
DataHolder::~DataHolder()
{
    delete [] len;
    delete [] ofs;
//    delete [] seq;
}

void DataHolder::get_len(const char *filename, int size)
{
    ifstream fin(filename);
    char s[1000];
    seq_num = 0;
    vector<bit64_t> len_vec;
    bit64_t cur = 0;

    while(fin.getline(s, 1000))
    {
        if(s[0] == '>')
        {
            if(size != -1 && seq_num == size) break;
            seq_num++;
            if(cur) len_vec.push_back(cur);
            cur = 0;
            continue;
        }
        cur += strlen(s);
    }
    if(cur) len_vec.push_back(cur), seq_num++;

//    ofs = new bit64_t[1000];
//    len = new bit64_t[1000];
    ofs = (bit64_t*) malloc((seq_num+1) * sizeof(bit64_t));
    len = (bit64_t*) malloc(seq_num * sizeof(bit64_t));

    ofs[0] = 0;
    for(int i = 0; i < seq_num; i++)
    {
        len[i] = len_vec[i];
        ofs[i+1] = ofs[i] + len_vec[i];
    }
    seq = new char[ofs[seq_num]];
//    if(seq_num == 193) seq = _seq; else seq = new char[ofs[seq_num]];
    fin.close();
}

void DataHolder::load(const char *filename, int size)
{
    get_len(filename, size);
    ifstream fin(filename);
    char s[1000];

    int seq_now = -1;
    bit64_t cur_ofs = 0;
    while(fin.getline(s, 1000))
    {
        if(s[0] == '>')
        {
            if(seq_now + 1 == seq_num) break;
            seq_now++;
            continue;
        }

        int tmp = strlen(s);
        strcpy(seq+cur_ofs, s);
        cur_ofs += tmp;
    }
    fin.close();
}

void DataHolder::get2(char *res, ref_id_t id, ref_loc_t start, ref_loc_t end)
{
//    cout << "in get2: " << id << " " << start << " " << end << endl;
    for(int i = start; i < end; i++) res[i-start] = seq[ofs[id] + i];
    res[end-start] = 0;
}

void DataHolder::get(char *res, ref_id_t id, ref_loc_t loc, bool if_reversed)
{
//    cout << "in get: " << id << " " << loc << " " << if_reversed << endl;
    if(if_reversed)
    {
        loc = len[id] - 1 - loc;
        if(loc < 2*_seed_size) get2(res, id, 0, _seed_size * 3);
        else if(loc > len[id]-_seed_size) get2(res, id, len[id]-3*_seed_size, len[id]);
        else get2(res, id, loc-2*_seed_size, loc+_seed_size);
        comp_reverse(res);
    }
    else
    {
        if(loc < _seed_size) get2(res, id, 0, _seed_size * 3);
        else if(loc > len[id]-_seed_size*2) get2(res, id, len[id]-3*_seed_size, len[id]);
        else get2(res, id, loc-_seed_size, loc+2*_seed_size);
    }
}

int DataHolder::get_hollow(char *res, ref_id_t id, ref_loc_t loc, bool if_reversed)
{
//    cout << "in get: " << id << " " << loc << " " << if_reversed << endl;
    if(if_reversed)
    {
        loc = len[id] - 1 - loc;
        int lower = max(0, (int)loc - 2 * _seed_size + 1);
        int upper = min(len[id], (bit64_t)loc + _seed_size + 1);
        get2(res, id, lower, upper);
        comp_reverse(res);
    }
    else
    {
        int lower = max(0, (int)loc - _seed_size);
        int upper = min(len[id], (bit64_t)loc + 2 * _seed_size);
        get2(res, id, lower, upper);
    }

    return loc < _seed_size ? loc : _seed_size;
}

// --------------------------------- StatConverter ---------------------------------
StatConverter::StatConverter(const char *filename)
{
    fin = ifstream(filename, ios_base::in | ios_base::binary);
}

StatConverter::~StatConverter()
{
    fin.close();
}

bool StatConverter::next_pair(bit64_t &read_val, bit64_t &ref_val)
{
    if(!fin.read((char *) &read_val, sizeof(read_val))) return false;
    if(!fin.read((char *) &ref_val, sizeof(read_val))) return false;
    return true;
}

void StatConverter::convert(DataHolder &ref, DataHolder &reads, const char *output_fn)
{
    ofstream fout(output_fn);
    bit64_t read_val, ref_val;
    char *ref_seq = new char[ref._seed_size*3 + 10];
    char *reads_seq = new char[reads._seed_size*3 + 10];
    bool if_reversed;
    ref_id_t id;
    ref_loc_t loc;

    while(next_pair(read_val, ref_val))
    {
//        cout << read_val << "\t" << ref_val << endl;
        inv_transform(read_val, id, loc, if_reversed);
        reads.get(reads_seq, id, loc, if_reversed);
        inv_transform(ref_val, id, loc, if_reversed);
        ref.get(ref_seq, id, loc, if_reversed);

//        cout << reads_seq << " " << ref_seq << endl;
        fout << reads_seq << " " << ref_seq << endl;

    }
    fout.close();
}

void StatConverter::convert_hollow(DataHolder &ref, DataHolder &reads, const char *output_fn)
{
    ofstream fout(output_fn);
    bit64_t read_val, ref_val;
    char *ref_seq = new char[ref._seed_size*3 + 10];
    char *reads_seq = new char[reads._seed_size*3 + 10];
    bool if_reversed;
    ref_id_t id;
    ref_loc_t loc;

    while(next_pair(read_val, ref_val))
    {
//        cout << read_val << "\t" << ref_val << endl;
        inv_transform(read_val, id, loc, if_reversed);
        int ofs1 = reads.get_hollow(reads_seq, id, loc, if_reversed);
        inv_transform(ref_val, id, loc, if_reversed);
        int ofs2 = ref.get_hollow(ref_seq, id, loc, if_reversed);

//        cout << reads_seq << " " << ref_seq << endl;
        fout << ofs1 << " " << ofs2 << endl;
        fout << reads_seq << " " << ref_seq << endl;
    }
    fout.close();
}