//--------------------------------------------------------------------------------------------------------
// insertion, deletion, mutation以外にduplicationとcontractionを考慮した編集距離(ed)の計算
// Reference: Tamar Pinhas, Shay Zakov, Dekel Tsur and Michal Ziv-Ukelson
// "Efficient edit distance with duplications and contractions” Algorithms for Molecular Biology, 8:27 (2013)
// To compile, perform: g++ -std=c++20 -Wall --pedantic-errors -o EDDC EDDC.cpp
//--------------------------------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <vector>
#include <climits>
using namespace std;

int ins(char a) { return 3; }
int del(char a) { return 3; }
int dup(char a) { return 2; }
int cont(char a) { return 2; }
// 塩基置換はKimuraの2-parameterモデルを使用
int alpha = 1;
int beta = 3;
int mut(char a, char b)
{
    if (a == b) return 0;
    else if ((a == 'A' && b == 'G') || (a == 'G' && b  == 'A') || 
            (a == 'C' && b == 'T') || (a == 'T' && b == 'C')) return alpha;
    else return beta;
}

struct EDDC
{
    public:
        EDDC(const string & s, const string & t)
            : s_(s), t_(t)
            {}

        int compute_edit_distance()
        {
            // DPテーブルのサイズを決める
            int len_s = s_.size();
            int len_t = t_.size();
            int num_alphabet = alphabet_.size();

            ed_s_to_empty_.assign          (len_s + 1, vector<int>(len_s + 1, 0));
            ed_s_to_alphabet_.assign       (num_alphabet, vector<vector<int>>(len_s + 1, vector<int>(len_s + 1, 0)));
            ed_s_to_alphabet_nongen_.assign(num_alphabet, vector<vector<int>>(len_s + 1, vector<int>(len_s + 1, 0)));
            ed_empty_to_t_.assign          (len_t + 1, vector<int>(len_t + 1, 0));
            ed_alphabet_to_t_.assign       (num_alphabet, vector<vector<int>>(len_t + 1, vector<int>(len_t + 1, 0)));
            ed_alphabet_to_t_nonred_.assign(num_alphabet, vector<vector<int>>(len_t + 1, vector<int>(len_t + 1, 0)));
            edt_.assign                    (num_alphabet, vector<vector<int>>(len_s + 1, vector<int>(len_t + 1, 0)));
            ed_.assign                     (len_s + 1, vector<int>(len_t + 1, 0));

            // Stage 1: source文字列とtarget文字列のいずれかが空文字 or 1文字の場合の編集距離を計算
            // DPテーブルの初期化
            for (int i = 0; i < len_s; i++) ed_s_to_empty_[i][i + 1] = del(s_[i]);
            for (int i = 0; i < len_t; i++) ed_empty_to_t_[i][i + 1] = ins(t_[i]);
            for (int k = 0; k < num_alphabet; k++)
            {
                for (int i = 0; i < len_s; i++) ed_s_to_alphabet_[k][i][i + 1] = mut(alphabet_[k], s_[i]);
                for (int i = 0; i < len_t; i++) ed_alphabet_to_t_[k][i][i + 1] = mut(alphabet_[k], t_[i]);
            }
            
            for (int j = 2; j <= len_t; j++)
            {
                for (int i = j - 2; i >= 0; i--)
                {
                    // Equation 3: alphabet_[k]の1文字スタートかつ最初の操作がmutでない場合
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        vector<int> tmp1(j - i - 1, INT_MAX);
                        vector<int> tmp2(j - i - 1, INT_MAX);
                        vector<int> tmp3(j - i - 1, INT_MAX);
                        for (int h = i + 1; h < j; h++)
                        {
                            tmp1[h - i - 1] = ed_alphabet_to_t_[k][i][h] + ed_empty_to_t_[h][j];
                            tmp2[h - i - 1] = ed_empty_to_t_[i][h] + ed_alphabet_to_t_[k][h][j];
                            tmp3[h - i - 1] = dup(alphabet_[k]) + ed_alphabet_to_t_[k][i][h] + ed_alphabet_to_t_[k][h][j];
                        }
                        int tmp1_min = *min_element(tmp1.begin(), tmp1.end());
                        int tmp2_min = *min_element(tmp2.begin(), tmp2.end());
                        int tmp3_min = *min_element(tmp3.begin(), tmp3.end());
                        ed_alphabet_to_t_nonred_[k][i][j] = min({tmp1_min, tmp2_min, tmp3_min});
                    }

                    // Equation 2: alphabet_[k]の1文字スタートかつ最初の操作がalphabet_[l]へのmutの場合
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        vector<int> tmp;
                        for (int l = 0; l < num_alphabet; l++)
                        {
                            tmp.push_back(mut(alphabet_[k], alphabet_[l]) + ed_alphabet_to_t_nonred_[l][i][j]);
                        }
                        ed_alphabet_to_t_[k][i][j] = *min_element(tmp.begin(), tmp.end());
                    }

                    // Equation 1: 空文字スタートの場合
                    vector<int> tmp;
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        tmp.push_back(ins(alphabet_[k]) + ed_alphabet_to_t_[k][i][j]);
                    }   
                    ed_empty_to_t_[i][j] = *min_element(tmp.begin(), tmp.end());
                }
            }

            for (int j = 2; j <= len_s; j++)
            {
                for (int i = j - 2; i >= 0; i--)
                {
                    // Equation 6 : alphabet_[k]の1文字で終わりかつ最後の操作がmutでない場合
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        vector<int> tmp1(j - i - 1, INT_MAX);
                        vector<int> tmp2(j - i - 1, INT_MAX);
                        vector<int> tmp3(j - i - 1, INT_MAX);
                        for (int h = i + 1; h < j; h++)
                        {
                            tmp1[h - i - 1] = ed_s_to_alphabet_[k][i][h] + ed_s_to_empty_[h][j];
                            tmp2[h - i - 1] = ed_s_to_empty_[i][h] + ed_s_to_alphabet_[k][h][j];
                            tmp3[h - i - 1] = cont(alphabet_[k]) + ed_s_to_alphabet_[k][i][h] + ed_s_to_alphabet_[k][h][j];
                        }
                        int tmp1_min = *min_element(tmp1.begin(), tmp1.end());
                        int tmp2_min = *min_element(tmp2.begin(), tmp2.end());
                        int tmp3_min = *min_element(tmp3.begin(), tmp3.end());
                        ed_s_to_alphabet_nongen_[k][i][j] = min({tmp1_min, tmp2_min, tmp3_min});
                    }

                    // Equation 5: alphabet_[k]の1文字で終わりかつ最後の操作がalphabet_[l]からのmutの場合
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        vector<int> tmp;
                        for (int l = 0; l < num_alphabet; l++)
                        {
                            tmp.push_back(mut(alphabet_[l], alphabet_[k]) + ed_s_to_alphabet_nongen_[l][i][j]);
                        }
                        ed_s_to_alphabet_[k][i][j] = *min_element(tmp.begin(), tmp.end());
                    }

                    // Equation 4: 空文字で終わりの場合
                    vector<int> tmp;
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        tmp.push_back(del(alphabet_[k]) + ed_s_to_alphabet_[k][i][j]);
                    }   
                    ed_s_to_empty_[i][j] = *min_element(tmp.begin(), tmp.end());
                }
            }

            // Stage 2: source文字列とtarget文字列のどちらも2文字以上の場合の編集距離を計算
            // s_[0]とt_[0]のalphabet_のインデックスを取得
            int s0_idx = 0;
            int t0_idx = 0;
            for (int i = 0; i < num_alphabet; i++)
            {
                if (s_[0] == alphabet_[i]) break;
                else s0_idx++;
            }
            for (int i = 0; i < num_alphabet; i++)
            {
                if (t_[0] == alphabet_[i]) break;
                else t0_idx++;
            }

            // DPテーブルの初期化
            ed_[0][0] = 0;
            for (int i = 1; i <= len_t; i++)
            {
                ed_[0][i] = ed_empty_to_t_[0][i];
                ed_[1][i] = ed_alphabet_to_t_[s0_idx][0][i];
            }
            for (int i = 1; i <= len_s; i++)
            {
                ed_[i][0] = ed_s_to_empty_[0][i];
                ed_[i][1] = ed_s_to_alphabet_[t0_idx][0][i];
            }
            // 論文には書いてないけどedt_の1行目もEquation 9で初期化しておく必要がある
            for (int j = 2; j <= len_t; j++)
            {
                for (int k = 0; k < num_alphabet; k++)
                {
                    vector<int> tmp(j - 1, INT_MAX);
                    for (int h = 1; h < j; h++)
                    {
                        tmp[h - 1] = ed_[1][h] + ed_alphabet_to_t_[k][h][j];
                    }
                    edt_[k][1][j] = *min_element(tmp.begin(), tmp.end());
                }
            }

            for (int j = 2; j <= len_t; j++)
            {
                for (int i = 2; i <= len_s; i++)
                {
                    // Equation 9: t_[0,j)の末尾だけalphabet1文字から生成されるようなs_[0,i)とt_[0,j)の編集パス
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        vector<int> tmp(j - 1, INT_MAX);
                        for (int h = 1; h < j; h++)
                        {
                            tmp[h - 1] = ed_[i][h] + ed_alphabet_to_t_[k][h][j]; // s_[0,i)がt_[0,h)に変換され、alphabet_[k]がt_[h,j)に変換される場合の編集距離
                        }
                        edt_[k][i][j] = *min_element(tmp.begin(), tmp.end());
                    }

                    // Equation 8: s_[0,i)とt_[0,j)の編集距離
                    vector<int> tmp((i - 1) * num_alphabet, INT_MAX);
                    for (int k = 0; k < num_alphabet; k++)
                    {
                        for (int h = 1; h < i; h++)
                        {
                            int ed1 = ed_s_to_alphabet_[k][0][i] + ed_alphabet_to_t_[k][0][j]; // s_[0,i)をalphabet_[k]に変換し、それをさらにt_[0,j)に変換するときの編集距離
                            int ed2 = edt_[k][h][j] + ed_s_to_alphabet_[k][h][i];              // s_[h,i)がalphabet_[k]に変換され、それがt_[0,j)の末尾になるようなs_[0,i)とt_[0,j)の編集パス
                            tmp[(k * (i - 1)) + (h - 1)] = min({ed1, ed2});
                        }
                    }
                    ed_[i][j] = *min_element(tmp.begin(), tmp.end());
                }
            }
            
            return ed_[len_s][len_t];
        }

        vector<vector<int>> &         get_ed_s_to_empty()           { return ed_s_to_empty_;           }
        vector<vector<vector<int>>> & get_ed_s_to_alphabet()        { return ed_s_to_alphabet_;        }
        vector<vector<vector<int>>> & get_ed_s_to_alphabet_nongen() { return ed_s_to_alphabet_nongen_; }
        vector<vector<int>> &         get_ed_empty_to_t()           { return ed_empty_to_t_;           }
        vector<vector<vector<int>>> & get_ed_alphabet_to_t()        { return ed_alphabet_to_t_;        }
        vector<vector<vector<int>>> & get_ed_alphabet_to_t_nonred() { return ed_alphabet_to_t_nonred_; } 
        vector<vector<vector<int>>> & get_edt()                     { return edt_;                     }
        vector<vector<int>> &         get_ed()                      { return ed_;                      }

    private:
        const string                s_;                       // source文字列
        const string                t_;                       // target文字列
        vector<char>                alphabet_ = {'A', 'C', 'G', 'T'};
        vector<vector<int>>         ed_s_to_empty_;           // s_[i, j]から空文字列への編集距離
        vector<vector<vector<int>>> ed_s_to_alphabet_;        // s_[i, j]からalphabet_[k]への編集距離
        vector<vector<vector<int>>> ed_s_to_alphabet_nongen_; // s_[i, j]からalphabet_[k]へのnon-generatingな操作による編集距離
        vector<vector<int>>         ed_empty_to_t_;           // 空文字列からt_[i, j]への編集距離
        vector<vector<vector<int>>> ed_alphabet_to_t_;        // alphabet_[k]からt_[i, j]への編集距離
        vector<vector<vector<int>>> ed_alphabet_to_t_nonred_; // alphabet_[k]からt_[i, j]へのnon-reducingな操作による編集距離
        vector<vector<vector<int>>> edt_;                     // alphabet_[k]を経由したs_[0, i]からt_[0, j]への編集距離
        vector<vector<int>>         ed_;                      // s_[0, i]からt_[0, j]への編集距離

        void print_dp_tables()
        {
            cout << "ED: S to Empty:" << "\n";
            for (auto & row : ed_s_to_empty_)
            {
                for (const auto & val : row)
                {
                    cout << val << " ";
                }
                cout << "\n";
            }
            
            cout << "ED: S to Alphabet:" << "\n";
            for (int i = 0; i < ed_s_to_alphabet_.size(); i++)
            {
                cout << "Alphabet " << alphabet_[i] << ":\n";
                for (const auto & row : ed_s_to_alphabet_[i])
                {
                    for (const auto & val : row)
                    {
                        cout << val << " ";
                    }
                    cout << "\n";
                }
            }

            cout << "ED: S to Alphabet non-gen:" << "\n";
            for (int i = 0; i < ed_s_to_alphabet_nongen_.size(); i++)
            {
                cout << "Alphabet " << alphabet_[i] << ":\n";
                for (const auto & row : ed_s_to_alphabet_nongen_[i])
                {
                    for (const auto & val : row)
                    {
                        cout << val << " ";
                    }
                    cout << "\n";
                }
            }

            cout << "ED: Empty to T:" << "\n";
            for (const auto & row : ed_empty_to_t_)
            {
                for (const auto & val : row)
                {
                    cout << val << " ";
                }
                cout << "\n";
            }

            cout << "ED: Alphabet to T:" << "\n";
            for (int i = 0; i < ed_alphabet_to_t_.size(); i++)
            {
                cout << "Alphabet " << alphabet_[i] << ":\n";
                for (const auto & row : ed_alphabet_to_t_[i])
                {
                    for (const auto & val : row)
                    {
                        cout << val << " ";
                    }
                    cout << "\n";
                }
            }

            cout << "ED: Alphabet to T non-reducing:" << "\n";
            for (int i = 0; i < ed_alphabet_to_t_nonred_.size(); i++)
            {
                cout << "Alphabet " << alphabet_[i] << ":\n";
                for (const auto & row : ed_alphabet_to_t_nonred_[i])
                {
                    for (const auto & val : row)
                    {
                        cout << val << " ";
                    }
                    cout << "\n";
                }
            }

            cout << "EDT:" << "\n";
            for (int i = 0; i < edt_.size(); i++)
            {
                cout << "Alphabet " << alphabet_[i] << ":\n";
                for (const auto & row : edt_[i])
                {
                    for (const auto & val : row)
                    {
                        cout << val << " ";
                    }
                    cout << "\n";
                }
            }

            cout << "ED:" << "\n";
            for (const auto & row : ed_)
            {
                for (const auto & val : row)
                {
                    cout << val << " ";
                }
                cout << "\n";
            }
        }
};

int main()
{
    string s = "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT";
    string t = "ACGTACGTACGT";
    EDDC eddc(s, t);
    int distance = eddc.compute_edit_distance();
    cout << "Edit Distance: " << distance << endl;

    return 0;
}