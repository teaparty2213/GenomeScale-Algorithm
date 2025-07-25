//--------------------------------------------------------------------------------------------------------
// suffix array (sa)をdoubling法で構築するクラス
// Reference: N. Jesper Larsson and Kunihiko Sadakane. 
// "Faster suffix sorting” Theoretical Computer Science, 387 (2007): 258-272.
// To compile, perform: g++ -std=c++20 -Wall --pedantic-errors -o SALS SALS.cpp
// debugする際の注意として、sa_は絶対値を取ること。isa_はisa_[sa_[i]]の形で参照すること。
// sa_の要素は負の値を持つことがある。これは、ある要素のsuffix array中の位置が確定した(ソート済みグループである)ことを示す。
//--------------------------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <random>
using namespace std;

struct SaLs
{
    public:
        SaLs(const int len_seq)
            : seq_(""),  len_seq_(len_seq),    sa_(len_seq, 0),  isa_(len_seq, 0) 
            { init_rng(); gen_random_seq(); }
        SaLs(const string & seq)
            : seq_(seq), len_seq_(seq.size()), sa_(len_seq_, 0), isa_(len_seq_, 0) 
            { init_rng(); }
        SaLs(const string & seq, const vector<char> & alphabet)
            : seq_(seq), len_seq_(seq.size()), sa_(len_seq_, 0), isa_(len_seq_, 0), alphabet_(alphabet) 
            { init_rng(); }
        
        void build_suffix_array()
        {
            create_alphabet_map();
            init_sa_and_isa();
            num_order_ = 1;

            while (num_sorted_groups_ < len_seq_ && num_order_ <= len_seq_) 
            // 2番目の条件がないとsa_[len_seq_ - 1] = 0のときに無限ループ(例: seq="TGGGCCCCA$")
            {
                // h-orderのソート済みグループを見つける
                int left_idx = 0; // unsorted groupの左端
                while (left_idx < len_seq_)
                {
                    if (sa_[left_idx] < 0) { left_idx++; continue; }

                    int right_idx = left_idx;
                    while (right_idx < (len_seq_ - 1) 
                        && isa_[abs(sa_[right_idx + 1])] == isa_[abs(sa_[left_idx])]) 
                    {
                        right_idx++;
                    }

                    // ソート済みグループのisa_とsa_を更新
                    ternary_split_quick_sort(left_idx, right_idx);
                    update_isa_and_sa(left_idx, right_idx);
                    left_idx = right_idx + 1;
                }
                num_order_ *= 2;
            }
            
            // sa_の要素は全て負なので-1倍して元に戻す
            for (int i = 0; i < len_seq_; i++)
            {
                if (sa_[i] < 0) sa_[i] *= -1;
            }
            is_valid_sa();
        }
        
        string &         get_seq()               { return seq_; }
        int              get_seq_len()           { return len_seq_; }
        vector<int> &    get_sa()                { return sa_; }
        vector<int> &    get_isa()               { return isa_; }
        vector<char> &   get_alphabet()          { return alphabet_; }
        map<char, int> & get_alphabet_map()      { return alphabet_map_; }
        int              get_num_sorted_groups() { return num_sorted_groups_; }
        int              get_num_order()         { return num_order_; }

    private:
        string         seq_ {""};                             // suffix arrayを構築する対象文字列
        int            len_seq_ {0};                          // seq_の長さ
        vector<int>    sa_;                                   // suffix array
        vector<int>    isa_;                                  // inverse suffix array
        vector<char>   alphabet_ = {'$', 'A', 'C', 'G', 'T'}; // アルファベット(デフォルトはDNAの4塩基と'$')
        map<char, int> alphabet_map_;                         // アルファベットを整数に対応させるmap
        int            num_sorted_groups_ {0};                // ソート済みグループの数
        int            num_order_ {0};                        // h-order
        mt19937        rng_;                                  // 乱数生成器
        
        void init_rng()
        {
            random_device seed_gen;
            auto seed = seed_gen();
            rng_.seed(seed);
        }

        void gen_random_seq()
        {
            uniform_int_distribution<int> d(1, alphabet_.size() - 1);
            for (int i = 0; i < len_seq_; i++)
            {
                if (i < len_seq_ - 1) seq_.push_back(alphabet_[d(rng_)]);
                else seq_.push_back(alphabet_[0]); // 最後の1文字だけが'$'
            }
        }

        // アルファベットを整数に対応させるマップの作成
        void create_alphabet_map()
        {
            for (int i = 0; i < alphabet_.size(); i++) alphabet_map_[alphabet_[i]] = i;
        }

        void update_isa_and_sa(const int left_idx, const int right_idx)
        {
            // update中に更新したisa_を参照するとうまくいかない(例: "GGGGGATTTCTTTCTTCTCAACGGGTACC$")
            vector<int> tmp_isa = isa_;

            // 代表元を決定してtmp_isaを更新
            int rep_idx = right_idx; // グループの代表元
            for (int i = right_idx; i >= left_idx; i--)
            {
                if (num_order_ == 0) // 0-orderのときはisa_がセットされていないのでseq_を参照
                {
                    if (seq_[sa_[i]] != seq_[sa_[rep_idx]])
                    {
                        rep_idx = i; // 異なるグループが来たら代表元を更新
                    }
                }
                else
                {
                    if (isa_[abs(sa_[i]) + num_order_] != isa_[abs(sa_[rep_idx]) + num_order_]) 
                    {
                        rep_idx = i;
                    }
                }
                tmp_isa[abs(sa_[i])] = rep_idx;
            }
            
            // sa_[i]がsorted groupなら-1倍して目印とし、ソート済みグループの数を更新する
            for (int i = left_idx; i <= right_idx; i++)
            {
                bool is_sorted = false;
                if (i == left_idx)
                {
                    if (isa_[abs(sa_[i])] != isa_[abs(sa_[i + 1])]) is_sorted = true;
                }
                else if (i == right_idx)
                {
                    if (isa_[abs(sa_[i])] != isa_[abs(sa_[i - 1])]) is_sorted = true;
                }
                else
                {
                    if (isa_[abs(sa_[i])] != isa_[abs(sa_[i + 1])] && 
                        isa_[abs(sa_[i])] != isa_[abs(sa_[i - 1])]) is_sorted = true;
                }
                if (is_sorted)
                {
                    sa_[i] *= -1;
                    num_sorted_groups_++;
                }
            }

            // 更新したisaを反映
            isa_ = tmp_isa;
        }

        void init_sa_and_isa()
        {
            int num_alphabet = alphabet_.size();
            vector<int> cnt_alphabet(num_alphabet, 0);

            // 各文字の出現回数を数える
            for (int i = 0; i < len_seq_; i++) cnt_alphabet[alphabet_map_[seq_[i]]]++;
            
            // 累積和に変換
            for (int i = num_alphabet - 1; i >= 0; i--)
            {
                if (i == num_alphabet - 1) cnt_alphabet[i] = len_seq_ - cnt_alphabet[i];
                else                       cnt_alphabet[i] = cnt_alphabet[i + 1] - cnt_alphabet[i];
            }

            // sa_を初期化
            for (int i = 0; i < len_seq_; i++)
            {
                int idx = alphabet_map_[seq_[i]];
                sa_[cnt_alphabet[idx]] = i;
                cnt_alphabet[idx]++;
            }

            // isa_を初期化してsa_を更新
            update_isa_and_sa(0, len_seq_ - 1);
        }

        void ternary_split_quick_sort(const int left_idx, const int right_idx)
        {
            if (left_idx >= right_idx) return;
            uniform_int_distribution<> d(left_idx, right_idx);
            int pivot = isa_[abs(sa_[d(rng_)]) + num_order_];

            vector<int> small;
            vector<int> equal;
            vector<int> large;
            for (int i = left_idx; i <= right_idx; i++)
            {
                int idx = abs(sa_[i]) + num_order_;
                if      (isa_[idx] < pivot) small.push_back(sa_[i]);
                else if (isa_[idx] > pivot) large.push_back(sa_[i]);
                else                        equal.push_back(sa_[i]);
            }
            
            // sa_の更新
            for (int i = left_idx; i <= right_idx; i++)
            {
                if (i < left_idx + small.size()) 
                {
                    sa_[i] = small[i - left_idx];
                }
                else if (left_idx + small.size() <= i && i < left_idx + small.size() + equal.size())
                {
                    sa_[i] = equal[i - left_idx - small.size()];
                }
                else 
                {
                    sa_[i] = large[i - left_idx - small.size() - equal.size()];
                }
            }

            // smallとlargeのサイズが0でないときは再帰
            if (small.size() > 0) ternary_split_quick_sort(left_idx, left_idx + small.size() - 1);
            if (large.size() > 0) ternary_split_quick_sort(right_idx - large.size() + 1, right_idx);
        }

        void is_valid_sa()
        {
            for (int i = 0; i < len_seq_ - 1; i++)
            {
                string s1 = seq_.substr(sa_[i]);
                string s2 = seq_.substr(sa_[i + 1]);
                if (s1 > s2)
                {
                    cout << "Invalid case found" << "\n";
                    cout << "seq: " << seq_ << "\n";
                    cout << "SA: ";
                    for (auto & sa : sa_) cout << sa << " ";
                    cout << "\n";
                    cout << "position: " << i << "\n";
                    cout << "s1: " << s1 << "\n";
                    cout << "s2: " << s2 << "\n";
                    break;
                }
            }
        }
};

int main()
{
    vector<int> len_seq = {0, 10, 73, 100, 240, 777, 1000, 3511, 10000};
    for (int i = 0; i < len_seq.size(); i++)
    {
        cout << "Testing length: " << len_seq[i] << "\n";
        for (int j = 1; j <= 100; j++)
        {
            SaLs test(len_seq[i]);
            //cout << test.get_seq() << "\n";
            test.build_suffix_array();
        }
    }
    return 0;
}