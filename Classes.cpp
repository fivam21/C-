// Inner product
// *******************************************
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <fstream>
#include<cstdlib>

using namespace std;

//const int n = 2;
//typedef int matrix[2][2];
//
//void input_matrix(matrix matrix1);
//void matrix_mult(matrix matrix1, matrix matrix2, matrix matrix3);
//void output_matrix(matrix matrix1);
//
//int main(){
//    
//    matrix MatrixA = {0}, MatrixB = {0}, MultMatrix = {0};
//    
//    cout << "Enter first matrix values: "<<endl;
//    input_matrix(MatrixA);
//    
//    cout << "Enter second matrix values: "<<endl;
//    input_matrix(MatrixB);
//    
//    matrix_mult(MatrixA, MatrixB, MultMatrix);
//    
//    output_matrix(MultMatrix);
//    
//    return 0;
//}
//
//
//void input_matrix(matrix matrixA) {
//    cout << "Enter values for the 2x2 matrix:" << endl;
//    for (int i = 0; i < 2; i++) {
//        for (int j = 0; j < 2; j++) {
//            cout << "Enter value for matrix[" << i << "][" << j << "]: ";
//            cin >> matrixA[i][j];
//        }
//        cout << endl;
//    }
//}
//
//
//void matrix_mult(matrix a, matrix b, matrix c){
//    
//    for (int i = 0; i < n; ++i) {
//        for (int k = 0; k < n; ++k) {
//            c[i][k] = 0; // Initialize the result matrix element to 0
//            
//            for (int j = 0; j < n; ++j) {
//                c[i][k] += a[i][j] * b[j][k]; // Perform matrix multiplication
//            }
//        }
//    }
//}
//
//
//void output_matrix(matrix matrixA){
//    cout << "Matrix Multiplication leads to: " <<endl;
//    for (int i=0; i<n; i++){
//        for (int j=0; j < n; j++){
//            cout << matrixA[i][j] << "  ";
//        }
//        cout << endl;
//    }
//}


// ASIAN CALL WITH ARITHMETIC AVERAGING
// *******************************************
//#include <vector>
//
//const float pi = 3.14152;
//const double OneOverRootTwoPi = 1/(sqrt(2*M_PI));
//const int N = 252;
//const int T = 1;
//const float dt = 1.0 / N;
//
//double normal(double mean, double std);
//double NormalDensity(double x);
//float max(float num1, float num2);
//float average(vector<float> a);
//
//class AsianOption{
//private:
//    int _strike;
//    float _rate, _vol, _stock;
//    
//public:
//    AsianOption(){
//        _stock = 0;
//        _strike = 0;
//        _rate = 0;
//        _vol = 0;
//    }
//    
//    AsianOption(float stock, int strike, float rate, float vol){
//        _stock = stock;
//        _strike = strike;
//        _rate = rate;
//        _vol = vol;
//    }
//    
//    ~AsianOption(){}
//    
//    vector<float> EM_stock_sim(float stock, float rate, float vol){
//        
//        double dX=0, dS=0;
//        vector<float> stock_arr(N, 0);
//        
//        stock_arr[0] = stock;
//        
//        for (unsigned int i=1; i<N; i++){
//            dX = normal(0.0,1.0) * sqrt(dt);
//            dS = stock * ((rate*dt) + (dX * vol));
//            stock += dS;
//            stock_arr[i] = stock;
//        }
//        
//        return stock_arr;
//    }
//
//    
//    
//    float AsianPricing(float stock, int strike, float rate, float vol){
//        
//        float final_value = 0;
//        vector<float> option_prices(N, 0);
//        
//        for (unsigned int i=0; i<N; i++){
//            
//            vector<float> stock_prices(N, 0);
//            float average_price = 0;
//            
//            // Simulating stock
//            stock_prices = EM_stock_sim(stock, vol, rate);
//            average_price = average(stock_prices);
//            option_prices[i] = max(average_price - strike, 0);
//        }
//        
//        final_value = exp(-rate*T) * average(option_prices);
//        
//        return final_value;
//        }
//    
//};
//
//
//int main(){
//    
//    srand((unsigned)time(NULL));
//    
//    AsianOption pricer;
////    int stock=0, strike=0;
////    float rate=0, vol=0;
//    float final_value = 0;
//    
////    cout << "Input option parameters: " <<endl;
////    cin >> stock >> stike >> rate >> vol;
////    final_value = pricer.AsianPricing(stock, strike, rate, vol);
//    
//    final_value = pricer.AsianPricing(200, 100, 0.05, 0.5);
//    cout << final_value <<endl;
//}
//
//
//double NormalDensity(double x){
//    return OneOverRootTwoPi*exp(-x*x/2);
//    }
//
//
////Polar Marsaglia Transformation
//double normal(double mean, double std){
//    static int iset = 0;
//    static double gset;
//    double fac, r, v1, v2;
//
//    // create two normally-distributed numbers
//    if (iset == 0){
//        r = 0;
//        do{ //compute two possibles
//        v1 = 2.0 * rand() / RAND_MAX - 1.0;
//        v2 = 2.0 * rand() / RAND_MAX - 1.0;
//        // they define radius
//        r = v1 * v1 + v2 * v2;
//            } while (r >= 1.0 || r == 0.0);
//
//        // in unit circle? if not try again
//        fac = sqrt((-2 * log(r)) / r); // Box-Muller transform
//        gset = (v1 * fac);
//        iset = 1; // save one
//        v2 = v2 * fac * std + mean; // scale and return one
//        return v2;}
//
//    else{
//        iset = 0;
//        return (gset * std) + mean;} //scale and return the saved one
//}
//
//float max(float num1, float num2){
//    return (num1 > num2 ? num1 : num2);
//}
//
//float average(vector<float> a){
//    
//    float sum = 0;
//    for (unsigned int i=0; i<N; i++){
//        sum += a[i];
//    }
//    
//    return sum / N;
//}


//  OVERLOAD OPERATOR, DETERMINANT and Inversion
// *******************************************

//const int n = 10; // Assuming n is a constant
//
//class Vector {
//private:
//    int _matrix_size;
//    vector<float> _matrix;
//
//public:
//    Vector(){
//        _matrix_size = 0;
//    }
//
//    Vector(int matrix_size, vector<float> a){
//        _matrix_size = matrix_size;
//        _matrix = a;
//    }
//
//    ~Vector(){}
//
//    friend Vector operator+(const Vector& a, const Vector& b) {
//
//        vector<float> c(a._matrix_size, 0.0); // Initialize result vector with zeros
//
//        for (int i = 0; i < a._matrix_size; ++i) {
//            c[i] = a._matrix[i] + b._matrix[i];
//        }
//
//        return Vector(a._matrix_size, c);
//    }
//    
//    void print() const {
//        cout << "Vector content:" << endl;
//        for (int i = 0; i < _matrix_size; ++i) {
//            cout << _matrix[i] << " ";
//        }
//        cout << endl;
//    }
//    
//};
//
//int main() {
//    // Example usage of the Vector class and operator+
//    vector<float> vec1 = {1.0, 2.0, 3.0};
//    vector<float> vec2 = {4.0, 5.0, 6.0};
//
//    Vector v1(3, vec1);
//    Vector v2(3, vec2);
//
//    Vector sum = v1 + v2;
//    sum.print();
//    return 0;
//}


// ************************************
//const int n = 2;

//class MatrixManip{
//private:
//    float _matrixA[n][n];
//    float _matrixB[n][n];
//    
//    
//public:
//    MatrixManip(){
//        
//    }
//    
//    
//    MatrixManip(float matrixA[n][n], float matrixB[n][n]) {
//        // Copy the elements of matrixA and matrixB to _matrixA and _matrixB respectively
//        for (int i = 0; i < n; ++i) {
//            for (int j = 0; j < n; ++j) {
//                _matrixA[i][j] = matrixA[i][j];
//                _matrixB[i][j] = matrixB[i][j];
//            }
//        }
//    }
//    
//    ~MatrixManip(){} // destructor
//    
//    
//    float determinent(float matrix[n][n]){
//        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
//    }
//    
//    
//    void inversion(float matrix[n][n], float inv_matrix[n][n]){
//        // [a b][c d]
//        //        invert = 1/det [d, -b][-c a]
//        float det = 0;
//        det = determinent(matrix);
//        
//        inv_matrix[0][0] = matrix[1][1]*(1/det);
//        inv_matrix[1][1] = matrix[0][0]*(1/det);
//        inv_matrix[1][0] = matrix[0][1]*-1*(1/det);
//        inv_matrix[0][1] = matrix[1][0]*-1*(1/det);
//    }
//
//};
//
//int main(){
//    float matrixA[n][n] = {{4,2}, {2,4}}, inv_matrix[n][n] = {0};
//    MatrixManip matt;
// 
//    cout << matt.determinent(matrixA) <<endl;
//    matt.inversion(matrixA, inv_matrix);
//    
//    cout << "Inverse matrix:" << endl;
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            cout << inv_matrix[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    // Checking inversion correct should = I
//    cout << "Inversion check: " <<endl;
//    for (unsigned int i = 0; i < n; ++i){
//        for (unsigned int k = 0; k < n; ++k){
//            float temp = 0;
//            for (unsigned int j = 0; j < n; ++j){
//                temp += inv_matrix[i][j] * matrixA[j][k];
//            }
//            cout << temp <<endl;
//        }
//    }
//    
//    return 0;
//}






// REVERSE A STRING
// *******************************************
//#include <cstring>
//
//const int n = 100;
//char* reverse_Str(char string[n]);
//
//int main() {
//    
//    char string[n];
//    
//    cout<< "Enter string: "<<endl;
//    cin.getline(string, n);
//    
//    cout << reverse_Str(string) << endl;
//
//    return 0;
//}
//
//char* reverse_Str(char string[n]) {
//    char* reversed_string = new char[strlen(string) + 1]; // Allocate memory for the reversed string (+1 for null terminator)
//    long i = 0;
//    
//    for (int i = 0; i < strlen(string); i++) { // Fix loop initialization and condition
//        reversed_string[i] = string[strlen(string) - i - 1]; // Correct array indexing
//    }
//    
//    // Add whitespace if the length of the reversed string is less than n
//    if (strlen(reversed_string) < n) {
//        for (long i = strlen(reversed_string); i < n; i++) {
//            reversed_string[i] = ' '; // Append whitespace
//        }
//        reversed_string[n] = '\0'; // Add null terminator
//    }
//    return reversed_string;
//}


// Product of 2**n, AUC f(x) = x**2 via MC
// *************************************************

//int main(){
//    int n = 3;
//    float prod = 1;
//    
////    cout.setf(ios::scientific);
//    cout.precision(20);
//    for (unsigned int i = 1; i <= n; i++){
//        prod *= pow(2,i);
//    }
//    cout << prod<<endl;
//}


// ******* TRAPEZIUM RULE ****************
// trapz = h(0.5*f0 + sumfn-1 + 0.5*fn)
//int main(){
//    float range1 = 0, range2 = 1, disc = 10; //MC_samples = 100,
//    float f0 = pow(range1, 2), fn = pow(range2, 2), area=0, sum=0;
//    
//    float h = (range2 - range1) / disc;
//    
//    //    for (unsigned int i=0; i<MC_samples; i++){
//    for (unsigned int i=1; i<disc; i++){
//        float x = pow(range1 + (i * h), 2);
//        sum += x;
//    }
//    
//    area = h*(sum+(0.5*f0)+(0.5*fn));
//    cout << "Trapezoidal area approximation: " << area <<endl;
//    
//    return 0;
//}


// ******* SIMPSON RULE ****************
// simp = h/3*(f0 + f2n + 4*sum_n(f2i-1) + 2*sum_n-1(f2i))
//int main(){
//    float range1=0, range2=1, area=0, sum1 = 0, sum2 = 0;
//    float disc = 10000, h=0, f0 = 0, f2n = pow(range2, 2);
//    
//    h = (range2 - range1)/disc;
//    
//    // odds
//    for (unsigned int i=1; i<=disc; i+=2){
//        float x = pow(range1 + i*h, 2);
//        sum1 += x;
//    }
//    
//    // evens
//    for (unsigned int i=2; i<=disc; i+=2){
//        float x = pow(range1 + i*h, 2);
//        sum2 += x;
//    }
//    
//    area = h/3*(4*sum1 + 2*sum2 + f0 + f2n);
//    cout << "Trapezoidal area approximation: " << area <<endl;
//
//    return 0;
//}


// Overloading the >>  and << operators
//class Patient{
//private:
//    string _name;
//    string _NHSnum;
//    string _GPname;
//
//public:
//    Patient(){
//        _name = "";
//        _NHSnum = "";
//        _GPname = "";
//    }
//
//    Patient(string name, string NHSnum, string GPname){
//        _name = name;
//        _NHSnum = NHSnum;
//        _GPname = GPname;
//    }
//
//    ~Patient(){}
//
//    string getname(){
//        return _name;}
//
//    string getNHSnum(){
//        return _NHSnum;}
//
//    string getGPname(){
//        return _GPname;}
//
//    void setname(string name){
//        _name = name;}
//
//    void setNHSnum(string NHSnum){
//        _NHSnum = NHSnum;}
//
//    void setGPname(string GPname){
//        _GPname = GPname;}
//};
//
//ostream& operator<<(ostream& os, class Patient& patient){
//    os << "Name: " << patient.getname() << "\n"
//       << "NHS number: " << patient.getNHSnum() << "\n"
//       << "GP name: " << patient.getGPname() << endl;
//    return os;
//}
//
//istream& operator>>(istream& is, class Patient& patient){
//    string name, nhsnum, gpname;
//    cout << "What the name: " << endl;
//    getline(is, name);
//    cout << "What the NHS: " << endl;
//    getline(is, nhsnum);
//    cout << "What the GP: " << endl;
//    getline(is, gpname);
//
//    patient.setname(name);
//    patient.setNHSnum(nhsnum);
//    patient.setGPname(gpname);
//
//    return is;
//}
//
//int main(){
//
//    Patient p("peter", "1290", "montalto");
//    cout << p;
//
//    cin >> p;
//    cout << p;
//
//    return 0;
//}



// A function swap which interchanges two doubles. Give two different
// implementations of this function, such that the swapping of the two real numbers
// reflects back in the caller of the function swap. Both function implementations should use pointers

//void swapDoubles1(double* a, double* b);
//void swapDoubles2(double* a, double* b);
//    
//int main() {
//    double x = 3.5, y = 7.2;
//    cout << "Before swap: x = " << x << ", y = " << y << endl;
//    swapDoubles1(&x, &y);
//    cout << "After swap: x = " << x << ", y = " << y << endl;
//    swapDoubles2(&x, &y);
//    cout << "After swap: x = " << x << ", y = " << y << endl;
//}
//
//void swapDoubles1(double* a, double* b) {
//    double temp = *a;
//    *a = *b;
//    *b = temp;
//}
//
//void swapDoubles2(double* a, double* b){
//    *a += *b;
//    *b = *a - *b;
//    *a -= *b;
//}



// Function
//void pos_av(vector<float> x, int n);
//
//int main(){
//    
//    vector<float> x;
//    unsigned int n;
//    
//    cout << "Enter Length of Array: " <<endl;
//    cin >> n;
//    x.resize(n);
//    cout << "Enter Array Values: " <<endl;
//    for (unsigned int i=1; i<=n; i++){
//        cin >> x[i-1];
//    }
//    
//    pos_av(x, n);
//    
//    return 0;
//}
//
//void pos_av(vector<float> x, int n){
//    float sum = 0, count = 0;
//    
//    for (unsigned int i=1; i<=n; i++){
//        if (x[i-1] > 0){
//            sum += x[i-1];
//            count += 1;
//        }
//    }
//    cout << sum/count<<endl;
//}


// A function sum that which takes a pointer to a
//function int funct(int) and return sum_N(funct(i))
//int square(int x) {
//    return x * x;
//}
//
//// Sum function taking a function pointer and an integer N
//int sum(int (*funct)(int), int N) {
//    int totalSum = 0;
//    for (int i = 1; i <= N; ++i) {
//        totalSum += funct(i); // Call the function pointed to by funct for each i
//    }
//    return totalSum;
//}
//
//int main() {
//    int N = 5; // Example: Calculate sum from 1 to 5 of the square of i
//
//    // Call the sum function with a pointer to the square function and integer N
//    int result = sum(square, N);
//    cout << "The sum of squares from 1 to " << N << " is " << result << endl;
//
//    return 0;
//}



////A function reverse which takes a pointer to a list of floats
////and the length of the list and reverses it. For example the list
//// (10, 9, 8, 7) should be changed to (7, 8, 9, 10) for example.
//#include <cstring>
//void reverse(float* a, int n);
//
//int main() {
//    float num[] = {10, 9, 8, 7};
//    int n = 4;
//    float* pointerA = num;
//    
//    // Call the reverse function
//    reverse(pointerA, n);
//    
//    // Output the reversed array
//    cout << "Reversed array: ";
//    for (int i = 0; i < n; i++) {
//        cout << num[i] << " ";
//    }
//    cout << endl;
//    
//    return 0;
//}
////
//void reverse(float* a, int n) {
//    float* start = a;
//    float* end = a + n - 1;
//
//    // Swap values while start is less than end
//    while (start < end) {
//        // Swap the values pointed to by start and end
//        float temp = *start;
//        *start = *end;
//        *end = temp;
//
//        // Move pointers towards each other
//        ++start;
//        --end;
//    }
//}


// DISCOUNTING
//float dis_compounding(float rate1, float rate2, float value, float years1, float years2);
//float cont_compounding(float rate1, float rate2, float value, float years1, float years2);
//
//int main(){
//    
//    cout<<"3yrs disc: "<<dis_compounding(0.05, 0, 100, 3, 0)<<endl;
//    cout<<"1.5yrs disc: "<<dis_compounding(0.04, 0.062, 100, 0.5, 1)<<endl;
//    cout<<"3yrs cont: "<<cont_compounding(0.05, 0, 100, 3, 0)<<endl;
//    cout<<"1.5yrs cont: "<<cont_compounding(0.04, 0.062, 100, 0.5, 1)<<endl;
//    return 0;
//}
//
//float dis_compounding(float rate1, float rate2, float value, float years1, float years2){
//    if (rate2 == 0){
//        return value / pow((1+rate1/years1), years1);
//    }
//    else{
//        float value1 = value / pow((1+rate2), years2);
//        return value1 / pow((1+(rate1)), years1);
//    }
//}
//
//float cont_compounding(float rate1, float rate2, float value, float years1, float years2){
//    if (rate2 == 0){
//        return value / exp(rate1*years1);
//    }
//    else{
//        float value1 = value / exp(rate2 * years2);
//        return value1 / exp(rate1 * years1);
//    }
//}



//Monte carlo in array for EM method for stock prices - Assume you have normal distribution function
//double normal(double, double); // function prototype
//
//int main()
//{
//    srand((unsigned)time(NULL));
//    ofstream print;
//    
//    print.open("results.xls");
//    
//    long N = 1000;
//    double asset = 100, IR = 0.05, vol = 0.2; // variables & parameters
//    double dt = 1.0 / N; // step-size for time
//    print << 0 << '\t' << asset << endl;
//    
//    for (unsigned short int i = 1; i <= N; i++) {
//        double time = i * dt;
//        double dX = normal(1.0, 0.0) * sqrt(dt);
//        double dS = asset * (IR * dt + vol * dX);
//        asset += dS;
//        print << time << '\t' << asset << endl;
//    }
//    
//    print.close();
//    
//    return 0;
//}