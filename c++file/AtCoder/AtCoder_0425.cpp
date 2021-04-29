#include <cstdlib>
#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#define rep(i, n) for (int i = 0; i < (int)(n); i++)

using namespace std;



void sort_2(vector<int> &v, int x){
    int n = v.size();
    int max=n,min=0;
    if(n==0 || v[n-1]<=x){
        v.push_back(x);
    }
    else if(v[0]>x){
        v.insert(v.begin(),x);
    }
    else{
        rep(i,10){
            if(v[(max+min)/2]<x){
                min = (max+min)/2;
            }
            else{
                max = (max+min)/2 + 1;
            }
            int ss = (int)pow(2.0,double(i));
            if(ss>=n){
                break;
            }
        }
        for (int i = min; i < max; i++){
            if(v[i]<=x && x<=v[i+1]){
                v.insert(v.begin()+i+1,x);
                break;
            }
        } 
    }
    /*
    rep(i,v.size()){
        cout << v[i];
    }
    cout << "\t" << min << "\t" << max << endl;*/
}

void SoD(){
    long long ans=0;
    int N;
    cin >> N;
    vector<int> v;
    rep(i,N){
        int x;
        cin >> x;
        sort_2(v,x);
    }
    rep(i,N){
        //cout << v[i];
        ans += (long long)(2*i -N +1) * v[i];
    }
    //cout << endl;
    cout << ans;
}

void SoD2(){
    long long ans=0;
    int N;
    cin >> N;
    vector<int> v(N);
    rep(i,N){
        cin >> v[i];
    }
    sort(v.begin(),v.end());
    rep(i,N){
        ans += (long long)(2*i -N +1) * v[i];
    }
    cout << ans;

}

int main(){
    SoD2();
    return 0;
}