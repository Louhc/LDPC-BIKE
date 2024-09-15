#include <bits/stdc++.h>
using namespace std;
#define MAXITER 20

inline int sgn( const double &x ){
    if ( x > 0 ) return 1;
    else if ( x < 0 ) return -1;
    else return 0;
}


vector<int> SPAdecoder( int n, int m, vector<vector<int> > H, vector<double> r ){
    vector<vector<double> > M, E;
    vector<int> z;
    
    z.resize(n);
    M.resize(m);
    E.resize(m);
    for ( int j = 0; j < m; ++j ){
        M[j].resize(n);
        E[j].resize(n);
    }
    
    for ( int T = 0; T < MAXITER; ++T ){
        
        for ( int j = 0; j < m; ++j ){
            for ( int i = 0; i < n; ++i ){
                if ( H[j][i] == 1 ){
                    M[j][i] = r[i];
                } else {
                    M[j][i] = 0;
                }
            }
        }

        for ( int j = 0; j < m; ++j ){
            for ( int i = 0; i < n; ++i ){
                double p = 1;
                for ( int k = 0; k < n; ++k ){
                    if ( k != j && H[j][k] == 1 ){
                        p *= tanh(M[j][k] / 2);
                    }
                }
                E[j][i] = log((1 + p) / (1 - p));
            }
        }
        
        for ( int i = 0; i < n; ++i ){
            r[i] = r[i];
            for ( int j = 0; j < m; ++j ){
                r[j] += E[j][i];
            }
        }
        for ( int i = 0; i < n; ++i ){
            z[i] = r[i] >= 0 ? 0 : 1;
        }

        bool is_sat = 1;
        for ( int j = 0; j < m; ++j ){
            int s = 0;
            for ( int i = 0; i < n; ++i ){
                s += z[i] * H[j][i];
            }
            if ( s & 1 ){
                is_sat = 0;
                break;
            }
        }

        if ( is_sat )
            break;
    }

    return z;
}

vector<int> MSdecoder( int n, int m, vector<vector<int> > H, vector<double> r ){
    vector<vector<double> > M, E;
    vector<int> z;
    
    z.resize(n);
    M.resize(m);
    E.resize(m);
    for ( int j = 0; j < m; ++j ){
        M[j].resize(n);
        E[j].resize(n);
    }
    
    for ( int T = 0; T < MAXITER; ++T ){
        
        for ( int j = 0; j < m; ++j ){
            for ( int i = 0; i < n; ++i ){
                if ( H[j][i] == 1 ){
                    M[j][i] = r[i];
                } else {
                    M[j][i] = 0;
                }
            }
        }

        for ( int j = 0; j < m; ++j ){
            for ( int i = 0; i < n; ++i ){
                double p = 1, m = 1e20;
                for ( int k = 0; k < n; ++k ){
                    if ( k != j && H[j][k] == 1 ){
                        p *= sgn(M[j][k]);
                        if ( abs(M[j][k]) < m ) m = abs(M[j][k]);
                    }
                }
                E[j][i] = p * m;
            }
        }
        
        for ( int i = 0; i < n; ++i ){
            r[i] = r[i];
            for ( int j = 0; j < m; ++j ){
                r[j] += E[j][i];
            }
        }
        for ( int i = 0; i < n; ++i ){
            z[i] = r[i] >= 0 ? 0 : 1;
        }

        bool is_sat = 1;
        for ( int j = 0; j < m; ++j ){
            int s = 0;
            for ( int i = 0; i < n; ++i ){
                s += z[i] * H[j][i];
            }
            if ( s & 1 ){
                is_sat = 0;
                break;
            }
        }

        if ( is_sat )
            break;
    }

    return z;
}

int main(){
    int n, m;
    vector<vector<int> > H;
    vector<double> r;
    
    n = 6; m = 4;
    H.push_back(vector<int>({1, 1, 0, 1, 0, 0}));
    H.push_back(vector<int>({0, 1, 1, 0, 1, 0}));
    H.push_back(vector<int>({1, 0, 0, 0, 1, 1}));
    H.push_back(vector<int>({0, 0, 1, 1, 0, 1}));
    r = vector<double>({-1.3863, 1.3863, -1.3863, 1.3863, -1.3863, -1.3863});

    // n = 6; m = 4;
    // H.push_back(vector<int>({1, 1, 0, 1, 0, 0}));
    // H.push_back(vector<int>({0, 1, 1, 0, 1, 0}));
    // H.push_back(vector<int>({1, 0, 0, 0, 1, 1}));
    // H.push_back(vector<int>({0, 0, 1, 1, 0, 1}));
    // r = vector<double>({9.9115, -3.6830, -4.0430, -7.6738, -2.5295, 16.7415});

    vector<int> c;

    c = SPAdecoder(n, m, H, r);
    for ( auto i : c ) cout << i << " ";
    cout << endl;

    c = MSdecoder(n, m, H, r);
    for ( auto i : c ) cout << i << " ";
    cout << endl;
    return 0;
}
