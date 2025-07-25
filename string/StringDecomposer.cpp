//---------------------------------------------------------------------------------------------------------------------------------
// 反復配列をユニットに分解する動的計画法アルゴリズム
// Reference: Tatiana Dvorkina, Andrey V. Bzikadze and Pavel A. Pevzner
// "The string decomposition problem and its applications to centromere analysis and assembly” Bioinformatics, 36 (2020): i93-i101.
// To compile, perform: g++ -std=c++20 -Wall --pedantic-errors -o StringDecomposer StringDecomposer.cpp
//---------------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <climits>
using namespace std;

const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -1;

int score(char a, char b)
{
    return (a == b) ? MATCH : MISMATCH;
}

struct StringDecomposer
{
    public:
        StringDecomposer(const string & seq, const vector<string> & blocks)
            : seq_(seq), blocks_(blocks), dp_(blocks.size())
            {}
        
        void decompose()
        {
            int len_seq = seq_.size();
            int num_blocks = blocks_.size();

            // dp[b][i][j]: blocks[b][0..i]とseq[0..j]の最適スコア
            for (int b = 0; b < num_blocks; ++b)
            {
                int len_block = blocks_[b].size();
                dp_[b].assign(len_block + 1, vector<int>(len_seq + 1, INT_MIN));
            }

            // dpテーブルの初期化
            for (int b = 0; b < num_blocks; ++b)
            {
                for (int i = 0; i <= blocks_[b].size(); ++i) dp_[b][i][0] = GAP * i;
                for (int j = 1; j <= len_seq;           ++j) dp_[b][0][j] = 0; // 各blockのページをgluedした部分
            }

            // dpテーブルを埋める
            for (int j = 1; j <= len_seq; ++j)
            {
                for (int b = 0; b < num_blocks; ++b)
                {
                    int len_block = blocks_[b].size();
                    for (int i = 1; i <= len_block; ++i)
                    {
                        int s_match = dp_[b][i - 1][j - 1] + score(blocks_[b][i - 1], seq_[j - 1]);
                        int s_del = dp_[b][i - 1][j] + GAP;
                        int s_ins = dp_[b][i][j - 1] + GAP;

                        dp_[b][i][j] = max({s_match, s_del, s_ins});
                    }
                }

                // block-switching edgeを考慮する
                vector<int> block_end_scores(num_blocks, 0);
                for (int b = 0; b < num_blocks; ++b)
                {
                    int len_block = blocks_[b].size();
                    block_end_scores[b] = dp_[b][len_block][j];
                }
                int max_end_score = *max_element(block_end_scores.begin(), block_end_scores.end());
                for (int b = 0; b < num_blocks; ++b) dp_[b][0][j] = max_end_score;
            }

            // sinkを求める
            int prev_score = INT_MIN;
            int sink_b = -1;
            int len_sink_block = -1;
            for (int b = 0; b < num_blocks; ++b)
            {
                int len_block = blocks_[b].size();
                prev_score = max(prev_score, dp_[b][len_block][len_seq]);
                if (prev_score == dp_[b][len_block][len_seq])
                {
                    sink_b = b;
                    len_sink_block = len_block;
                }
            }

            // trace backで最適パスを求める
            int b = sink_b;
            int i = len_sink_block;
            int j = len_seq;
            path_.emplace_back(b, i, j);
            while (i != 0 || j != 0)
            {
                if (prev_score == dp_[b][i - 1][j - 1] + score(blocks_[b][i - 1], seq_[j - 1]))
                { 
                    prev_score = dp_[b][i - 1][j - 1];
                    --i; --j;
                }
                else if (prev_score == dp_[b][i - 1][j] + GAP)
                {
                    prev_score = dp_[b][i - 1][j];
                    --i;
                }
                else if (prev_score == dp_[b][i][j - 1] + GAP)
                {
                    prev_score = dp_[b][i][j - 1];
                    --j;
                }
                path_.emplace_back(b, i, j);

                // glued部分に到達したらblock-switching edgeを遡る
                if (i == 0 && j > 0)
                {
                    for (int prev_b = 0; prev_b < num_blocks; ++prev_b)
                    {
                        int prev_len_block = blocks_[prev_b].size();
                        if (dp_[prev_b][prev_len_block][j] == prev_score)
                        {
                            b = prev_b;
                            i = blocks_[b].size();
                            path_.emplace_back(b, i, j);
                            break;
                        }
                    }
                }
            }
            reverse(path_.begin(), path_.end());

            // seq_の最適な分解を求める
            int block_start_idx = 0;
            for (int i = 1; i < path_.size(); ++i)
            {
                auto [prev_b, prev_i, prev_j] = path_[i - 1];
                auto [curr_b, curr_i, curr_j] = path_[i];

                if (curr_i == 0 || i == path_.size() - 1)
                {
                    int block_end_idx = max(prev_j, curr_j);
                    decomp_.push_back(seq_.substr(block_start_idx, block_end_idx - block_start_idx));
                    block_start_idx = block_end_idx;
                }
            }
        }

        vector<vector<vector<int>>>  & get_dp()     { return dp_;     }
        vector<tuple<int, int, int>> & get_path()   { return path_;   }
        vector<string>               & get_decomp() { return decomp_; }

    private:
        const string                  seq_;
        vector<string>                blocks_;
        vector<vector<vector<int>>>   dp_;         // dpテーブル
        vector<tuple<int, int, int>>  path_;       // 最適パスのblockインデックス, i, j
        vector<string>                decomp_;     // seq_の最適な分解
};

int main()
{
    /*string seq = "ACGTCGC";
    vector<string> blocks = {"ACGT", "ATAT", "CGCG"};*/
    string seq = "ACGTACGTACCTACGTTCGTACGT";
    vector<string> blocks = {"ACGT", "ACCT", "ACGTT", "TCGT"};

    StringDecomposer test(seq, blocks);
    test.decompose();

    /*vector<vector<vector<int>>>  & dp = test.get_dp();
    for (int b = 0; b < blocks.size(); ++b)
    {
        cout << "Block " << b << ":\n";
        for (int i = 0; i <= blocks[b].size(); ++i)
        {
            for (int j = 0; j <= seq.size(); ++j)
            {
                cout << dp[b][i][j] << " ";
            }
            cout << "\n";
        }
        cout << "--------------------\n";
    }*/
    
    /*vector<tuple<int, int, int>> & path = test.get_path();
    for (const auto & p : path)
    {
        cout << "Block: " << get<0>(p) << ", i: " << get<1>(p) << ", j: " << get<2>(p) << "\n";
    }*/

    vector<string> & decomp = test.get_decomp();
    for (auto & d : decomp)
    {
        cout << d << " ";
    }
    cout << "\n";
    return 0;
}