#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>

using namespace std;

vector<vector<double>> readCSV(const string& filename, char delimiter = ',') {
    vector<vector<double>> result;

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return result;
    }

    string line;
    while (getline(file, line)) {
        vector<double> row;
        stringstream ss(line);
        string token;

        while (getline(ss, token, delimiter)) {
            try {
                double value = stod(token);
                row.push_back(value);
            } catch (const invalid_argument& e) {
                cerr << "Invalid argument: " << e.what() << endl;
            } catch (const out_of_range& e) {
                cerr << "Out of range: " << e.what() << endl;
            }
        }

        result.push_back(row);
    }

    file.close();
    return result;
}


vector<vector<double>> matrix_multipy(vector<vector<double>> &m1, vector<vector<double>> &m2)
{

    vector<vector<double>> result(m1.size(), vector<double>(m2[0].size(), 0));

    for (size_t i = 0; i < m1.size(); i++)
    {
        for (size_t j = 0; j < m2[0].size(); j++)
        {
            for (size_t k = 0; k < m2.size(); k++)
            {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return result;
}

vector<vector<double>> matrix_divide(vector<vector<double>> &m1, vector<vector<double>> &m2)
{

    vector<vector<double>> result(m1.size(), vector<double>(m2[0].size(), 0));

    for (size_t i = 0; i < m1.size(); i++)
    {
        for (size_t j = 0; j < m1[0].size(); j++)
        {
            result[i][j] = m1[i][j] / m2[i][j];
        }
    }

    return result;
}

vector<vector<double>> matrix_transpose(vector<vector<double>> &m1)
{

    vector<vector<double>> result(m1[0].size(), vector<double>(m1.size(), 0));

    for (size_t i = 0; i < m1.size(); i++)
    {
        for (size_t j = 0; j < m1[0].size(); j++)
        {
            result[j][i] = m1[i][j];
        }
    }

    return result;
}

vector<vector<double>> matrix_dot(vector<vector<double>> &m1, vector<vector<double>> &m2)
{

    vector<vector<double>> result(m1.size(), vector<double>(m2[0].size(), 0));

    for (size_t i = 0; i < m1.size(); i++)
    {
        for (size_t j = 0; j < m1[0].size(); j++)
        {
            result[i][j] = m1[i][j] * m2[i][j];
        }
    }

    return result;
}

vector<vector<double>> matrix_addition(vector<vector<double>> &m1, vector<vector<double>> &m2, bool add_mode)
{

    vector<vector<double>> result(m1.size(), vector<double>(m2[0].size(), 0));

    for (size_t i = 0; i < m1.size(); i++)
    {
        for (size_t j = 0; j < m1[0].size(); j++)
        {
            result[i][j] = m1[i][j] + m2[i][j] * (-1 + 2 * add_mode);
        }
    }

    return result;
}

vector<vector<double>> matrix_rnd_populate(int n, int m)
{

    random_device rnd_device;

    mt19937_64 gen(rnd_device());
    uniform_real_distribution<double> distribution(0.0, 1.0);

    vector<vector<double>> result(n, vector<double>(m, 0));

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            result[i][j] = abs(distribution(gen));
        }
    }

    return result;
}

void matrix_print(vector<vector<double>> &m1)
{
    cout << "\n";

    for (int i = 0; i < m1.size(); i++)
    {
        for (int j = 0; j < m1[i].size(); j++)
        {
            cout << m1[i][j];
            cout << " ";
        }
        cout << "\n";
    }
}

bool matrix_threshold(vector<vector<double>> &m1, double threshold, bool cummulative)
{
    bool valid = false;

    if (cummulative)
    {
        double c = 0;

        for (int i = 0; i < m1.size(); i++)
        {
            for (int j = 0; j < m1[i].size(); j++)
            {
                c += abs(m1[i][j]);
            }
        }
        cout << c << "\n";
        if (c <= threshold)
        {
            valid = true;
        }
    }
    else
    {
        valid = true;

        for (int i = 0; i < m1.size(); i++)
        {
            for (int j = 0; j < m1[i].size(); j++)
            {
                if (abs(m1[i][j]) > threshold)
                {
                    valid = false;
                }
            }
        }
    }

    return valid;
}

int main()
{
    /*
    vector<vector<double>> H {{1, 2, 3, 4, 5, 6}, {7, 8, 9, 10, 11, 12}};
    vector<vector<double>> W {{1, 2}, {3, 4}, {5, 6}, {7, 8}};

    

    int V_rows = 4;
    int V_cols = 6;
    */
   
    int N_alg = 2;

    
    //H = matrix_rnd_populate(N_alg,V_cols);
    //W = matrix_rnd_populate(V_rows,N_alg);

    //vector<vector<double>> V1 = matrix_multipy(W, H);
    

    string filename = "C:/Users/Pedro/Desktop/NNMF/X_TABLE_2.txt";  // Replace with your CSV file name
    char delimiter = ',';

    vector<vector<double>> V = readCSV(filename, delimiter);

    vector<vector<double>> H = readCSV("C:/Users/Pedro/Desktop/NNMF/h.txt", delimiter);;
    vector<vector<double>> W = readCSV("C:/Users/Pedro/Desktop/NNMF/w.txt", delimiter);;

    H = matrix_rnd_populate(N_alg,V[0].size());
    W = matrix_rnd_populate(V.size(),N_alg);

    //matrix_print(V);
    //matrix_print(V1);

    //matrix_print(H);
    //matrix_print(W);

    vector<vector<double>> W_trans;
    vector<vector<double>> W_trans_W;
    vector<vector<double>> H_mult_num;
    vector<vector<double>> H_mult_den;
    vector<vector<double>> H_mult;
    vector<vector<double>> H_trans;
    vector<vector<double>> W_H;
    vector<vector<double>> W_mult_num;
    vector<vector<double>> W_mult_den;
    vector<vector<double>> W_mult;
    vector<vector<double>> V_final;
    vector<vector<double>> diff;

    cout << "test";

    for (size_t i = 0; i < 100000; i++)
    {

        W_trans = matrix_transpose(W);
        W_trans_W = matrix_multipy(W_trans, W);

        H_mult_num = matrix_multipy(W_trans, V);
        H_mult_den = matrix_multipy(W_trans_W, H);
        H_mult = matrix_divide(H_mult_num, H_mult_den);
        H = matrix_dot(H, H_mult);

        H_trans = matrix_transpose(H);
        W_H = matrix_multipy(W, H);

        W_mult_num = matrix_multipy(V, H_trans);
        W_mult_den = matrix_multipy(W_H, H_trans);
        W_mult = matrix_divide(W_mult_num, W_mult_den);
        W = matrix_dot(W, W_mult);

        V_final = matrix_multipy(W, H);
        diff = matrix_addition(V, V_final, false);

        if(matrix_threshold(diff, 0.1, true))
        {
            cout << "GOOD w/ " << i+1 << " iterations.";
            break;
        }

    }

    //matrix_print(H_mult);
    //matrix_print(W_mult);
    //matrix_print(V);
    //matrix_print(V_final);
    //matrix_print(diff);

    return 0;
}