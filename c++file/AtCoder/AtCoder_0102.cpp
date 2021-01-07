#include <cstdlib>
#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#define rep(i, n) for (int i = 0; i < (int)(n); i++)

using namespace std;

void LD(){
    int A,B;
    cin >> A >> B;
    int a[3],b[3];
    int sum_A=0,sum_B=0;
    rep(i,3){
        int P = (int)pow((double)10,(double)(i+1));
        int P0 = (int)pow((double)10,(double)i);
        a[i] = A%P;
        A = A - a[i];
        a[i] = a[i]/P0;
        b[i] = B%P;
        B = B - b[i];
        b[i] = b[i]/P0;
        sum_A += a[i];
        sum_B += b[i];
    }
    if(sum_A>=sum_B){
        cout << sum_A << endl;
    }
    else{
        cout << sum_B << endl;
    }

}


void GP(){
    int N;
    cin >> N;
    int x[N],y[N];
    rep(i,N){
        cin >> x[i] >> y[i];
    }
    int ans=0;
    rep(i,N){
        rep(j,i){
            double P = abs(y[i]-y[j])/((double)abs(x[i]-x[j]));
            if(abs(P)<=1){
                ans++;
            }
        }
    }
    cout << ans << endl;
}

void SAT(){
    int N;
    cin >> N;
    vector<string> s,t;
    rep(i,N){
        string test;
        cin >> test;
        if(test[0]=='!'){
            test.erase(0,1);
            t.push_back(test);
        }
        else{
            s.push_back(test);
        }
    }
    sort(s.begin(),s.end());
    sort(t.begin(),t.end());

    int count = 0;
    
    if(s.size()==0){
        cout << "satisfiable" << endl;
    }
    else{
        rep(i,t.size()){
            int s_min=0,s_max=s.size();
            rep(j,15){
                if(s[(s_min+s_max)/2]>t[i]){
                    s_max = (s_min+s_max)/2 + 1;
                }
                else{
                    s_min = (s_min+s_max)/2;
                }
            }
            for (int j = s_min; j < s_max; j++){
                if(t[i]==s[j]){
                    cout << s[j] <<endl;
                    count++;
                    break;
                }
            }
            if(count > 0){
                break;
            }
            
        }

        if(count==0){
            cout << "satisfiable" << endl;
        }
    }
    
}

void CM(){
    int N;
    cin >> N;
    long A[N],T[N];
    vector<long> C;
    long long a_sum=0;
    rep(i,N){
        cin >> A[i] >> T[i];
        C.push_back(2*A[i] + T[i]);
        a_sum += A[i];
    }
    sort(C.begin(),C.end());
    int count=0;
    rep(i,N){
        int X = N-i-1;
        a_sum -= C[X];
        count++;
        if(a_sum <0){
            break;
        }
    }
    cout << count << endl;
}

int main(){
    CM();
    return 0;
}