// left rotation in array 
//  [2,4,1,9] => [1,9,2,4] , k = 2
#include<bits/stdc++.h>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);
    
    int n;
    cin >> n;
    int arr[n];
    for(int i =0;i<n;i++){
        cin >> arr[i];
    }
       int k;
       cin >> k;
       for(int i = k;i<k+n;i++){
           cout << arr[i%n] << " ";
       }
       return 0;
}
