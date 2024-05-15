#include <iostream>
#include <cmath>

using namespace std;


// b^2-4ac code:
//****************************
//int main(){
//    float a=0, b=0, c=0, root1, root2, x;
//    
//    cout<<"Enter values a,b,c: "<<endl;
//    cin>>a;
//    cin>>b;
//    cin>>c;
//    
//    x = sqrt(b*b - 4*a*c);
//    root1 = (-b+x)/2*a;
//    root2 = (-b-x)/2*a;
//    cout<<"root of equation is "<< root1 <<"and"<<root2 <<endl;
//    return 0;}



// accrued interest:
//****************************
//int main(){
//    float interest, tv, acc_int, notional;
//    int years;
//    cout<<"interest rate: "<<endl;
//    cin>>interest;
//    cout<<"many years: "<<endl;
//    cin>>years;
//    cout<<"amount invested: "<<endl;
//    cin>>notional;
//    
//    acc_int = (notional * pow(1+interest, years)) - notional;
//    cout<<acc_int<<endl;
//    cout<<sizeof(acc_int)<<endl;
//    return 0;}


// The time period T for a simple pendulum can be determined from 𝑇 = 2𝜋sqrt(l\g) , where l is the
// length of string and g is acceleration due to gravity (9.81 ms^-2). Write a program, which
// calculates values of T for varying lengths of string.
//****************************
//int main(){
//    float len_str = 0, g = 9.81;
//    const float pi = 3.14;
//    string string;
//    
//    cout<<"enter g: "<<endl;
//    cin>>string;
//    len_str = size(string);
//    
//    cout<<"T: " << 2*pi*sqrt(len_str/g) <<endl;
//    cout<< len_str<<endl;
//    cout<< pi<<endl;
//    cout<< g<<endl;
//    return 0;}


//****************************
// Write a program that calculates the current value of future cash flows. For example consider
// the following investments: 100k 1.5yrs at IR: 4% first 6mth, 6.2% after
// Use both discrete and continuous compounding in seprate programs.
//****************************
//int main(){
//    int notional = 100000;
//    float int_rate1 = 0.04, int_rate2 = 0.062, first_disc=0, first_cont=0;
//    
//    first_disc = notional*pow((1+int_rate1), 0.5); //V*(1+int)^time
//    first_cont =notional*exp(int_rate1*0.5); //V*exp(rt)
//    
//    cout<< "Disc Comp: " << first_disc*pow((1+int_rate2), 1) <<endl;
//    cout<< "Cont Comp: " << first_cont*exp(int_rate2*1) <<endl;
//    
//    return 0;}


// outputing files to excel
//****************************
#include <fstream> // library containing files
//
//int main(){
//    ofstream out; // create new object (out) belonging to class ofstream for printing to a file, this is the users choice
//    out.open("table.xls"); /* out will be used instead of cout to output results to an excel file called table.xls */
//    for (short unsigned int t = 0; t <= 10; t++ ){
//        out << t << '\t'<< t*t <<endl ;// every time out is used, results are stored in table.xls
//        }
//        out.close(); // file table.xls is now closed
//    return 0;}


//Compute sum i^2 = n(n+1)(n+2)/6 for both sides
//****************************
// Rewrite the code to sum odds and evens for each of the above.
//int main(){
//    int n=0,i_n=0, odd_count=0, even_count=0;
//    unsigned int i;
//    
//    cout<<"n is: "<<endl;
//    cin>>n;
//    
//    i_n = (n*(n+1)*(2*n+1))/6;
//    for (i=1; i<=n;++i){
//        if (i%2==0){
//            even_count += pow(i,2);
//        }
//        else{
//            odd_count += pow(i,2);
//        }
//    }
//    
//    cout<< odd_count << " plus "<< even_count << " should equal " << i_n <<endl;
//    return 0;
//}



// conditional operator 
//****************************
//int main(){
//    float output_=0, x=0;
//    cout<<"x"<<endl;
//    cin>>x;
//    output_ = (x>=0) ? sqrt(x): sqrt(-x);
//    cout<<output_<<endl;}



// Discretise x^2 over 0.1 steps and use 4dp precision
//****************************
//#include <fstream>
#include <iomanip>

//int main(){
//    double x=0, i=1;
//    cout.setf(ios::fixed);
//    cout.setf(ios::scientific);
//    cout.precision(4);
//    cout.width(10);
//    for (i=1; i<=2; i+=0.1){
//        cout<< pow(i,2)<<endl;}
//    return 0;}


//Arrays - calc mu and std of set of numbers (up to 10)
//****************************
//const int MAX_NUMS = 10;
//
//float Average(float arr[], int len_arr);
//float st_devi(float arr[], int len_arr, float mu);
//
//int main(){
//    float mean=0, std=0;
//    int n=0;
//    
//    cout<<"Whats len of array"<<endl;
//    cin>>n;
//    
//    float arr_one[MAX_NUMS] = {0};
//    
//    cout<<"enter elements: "<<endl;
//    for (int i=0;i<n;++i){
//        cin>>arr_one[i];}
//    cout<< "Mean: " << Average(arr_one, n) << "\n"<<endl;
//    cout<< "STD: "<< st_devi(arr_one, n, mean) <<endl;
//    return 0;}
//
//float Average(float arr[], int len_arr){
//    float sum=0;
//    
//    for (int i=0;i<len_arr;++i){
//        sum += arr[i];
//    }
//    return sum/len_arr;}
//
//float st_devi(float arr[], int len_arr, float mu){
//    float sum=0;
//    
//    for (int i=0; i<len_arr;++i){
//        sum += pow((arr[i] - mu) ,2);
//    }
//    return sqrt(sum)/(len_arr);}


// Strings
//****************************
#include <cstring>
//using namespace std;
//int main(){
//    char name[30] = "That is it "; /* Initialised a string with 30 slots, ie array holding characters to a message. */
//    char another_string[80]; /* Declaration of another to hold 80 chars, not initialised. */
//    
//    strcpy(name, "Riaz"); /* This takes string “name”, and copies another string over it including sentinel O*/
//    
//    cout << name; /* To print out string, do not need [] */
//    
//    cout << " is of length " << strlen(name) << endl; /* This built in func. counts length of string, ignoring \O */
//    
//    cout << "Enter a name: ";
//    cin.getline(another_string, 80); /* cin >> stops executing at a space, so cin acts on 
//                                      func getline which takes the name of string together with its max length*/
//    cout << another_string << endl;
//    if ((strcmp(name, another_string))) /* Compares 2 strings to check equivalence, if so returns False (0) */
//        cout << "Goodbye"; // comes here when different strings
//    else
//        cout << "hello" << endl; // comes here when same
//    cout << strcat(name, another_string); /* This attaches 2nd string on to the end of 1st, without gap */
//    return 0;
//}




// Pointers
//****************************
//typedef int * IntPtr;
//int main(){
//    IntPtr ptr_a, ptr_b;
//    int num_c = 4, num_d = 7;
//    
//    ptr_a = &num_c; /* LINE 11 */ // setting pointer address to num_c
//    ptr_b = ptr_a; /* LINE 12 */ // both pointers pointing to num_c
//    
//    cout << *ptr_a << " " << *ptr_b << "\n"; //4, 4 (* = variable the pointer is looking at)
//    
//    ptr_b = &num_d; /* LINE 16 */ //pointer reassigned to num_d
//    
//    cout << *ptr_a << " " << *ptr_b << "\n"; //4, 7
//    
//    *ptr_a = *ptr_b; /* LINE 20 */ //pointer a now pointing at 7 too
//    
//    cout << *ptr_a << " " << *ptr_b << "\n"; //7, 7
//    cout << num_c << " " << *&*&*&num_c << "\n"; // 7, 7
//    return 0;
//}

//void assign_new_int(IntPtrType& ptr){
//    ptr = new int;
//        if (ptr == NULL){
//        cout << "Sorry, ran out of memory";
//        exit(1);
//        }
//    }



// Pointer arthmetic
//****************************
//int main(){
//    // an array with 5 elements.
//    double T[5] = { 1.3, 6.9, 8.4, 16.3, 50.0 };
//    double* p;
//    p = T;
//    // output each array element's value
//    cout << "Array values using pointer " << endl;
//    for (int i = 0; i < 5; i++)
//        {
//        cout << "*(p + " << i << ") : ";
//        cout << *(p + i) << endl;
//        }
//    
//    cout << "Array values using T as address " << endl;
//    for (int i = 0; i < 5; i++)
//        {
//        cout << "*(T + " << i << ") : ";
//        cout << *(T + i) << endl;
//        }
//    
//    cout << T<<'\t'<<p<<endl; // both give the same base address
//    cout << *T<<'\t' << *p << endl; // both the same dereference (first element of array)
//    cout << *T << '\t'<<*T + 1 << endl; // basic arithmetic //1.3    2.3
//    return 0;}




// Passing by pointer
//****************************
//void add_to_byvalue(double x);
//void add_to_bypointer(double* x);
//
//int main() {
//    double z = 4.5;
//    add_to_byvalue(z); // pass z
//    cout << z << endl; // ’4.5’
//    add_to_bypointer(&z); // pass &z, ie the address
//    cout << z << endl; // ’5.5’
//    return 0;}
//
//void add_to_byvalue(double x){
//x++;} // does not change x outside
//
//void add_to_bypointer(double* x){
//    (*x)++;} // does change x outside MUST BE WRAPPED IN ()


// Whats the output? Pointer arthmetic
//****************************
//int main(){
//    typedef int* intptr;
//    intptr ptr_a;
//    int arr[]={10,11,12,13,14,15,16,17,18,19};
//    
//    ptr_a=&arr[3]; //address of 13 [index 4]
//    cout<<(ptr_a+4)<<endl;
//    
//    cout<<ptr_a<<'\t'<<*ptr_a<<endl;
//    
//    ++ptr_a;
//    cout<<ptr_a<<'\t'<<*ptr_a<<endl;
//    
//    ptr_a+=3;
//    cout<<ptr_a<<'\t'<<*ptr_a<<endl;
//    
//    ptr_a-=6;
//    cout<<ptr_a<<'\t'<<*ptr_a<<endl;
//    
//    --ptr_a;
//    cout<<ptr_a<<'\t'<<*ptr_a<<endl;
//    
//    for(int i=0; i<=9; ++i){
//        cout<<ptr_a[i]<<endl;}
//    return 0;}


// Recursion
//****************************
//int sum_nums(int n);
//int main(){
//    int sum, n;
//    cout<<"Sum of numbers 1 to n: what's n:"<<endl;
//    cin>>n;
//    sum = sum_nums(n);
//    cout<< "Sum is " << sum <<endl;
//    return 0;}
//
//int sum_nums(int n){
//    if (n>0){
//        n += sum_nums(n-1);
//    }
//    return n;}


// Templates
//****************************
//template <class T>
//T myMin(T a, T b){
//    return (a < b ? a : b);
//}
//
//int main(){
//    typedef unsigned short int ushort;
//    
//    int int1(-13), int2(-24); // brackets same as =
//    float x=2.3, y=3.6;
//    ushort a(2), b(3);
//    
//    cout<<myMin(int1, int2)<<endl;
//    cout << myMin(x, y) << endl;
//    cout << myMin(a, b) << endl;
//    
//    return 0;
//}


// Time class
//****************************
//typedef unsigned short int ushort;
//
//class Time{
//private:
//    int _hours, _mins, _secs;
//
//public:
//    Time(){ // default (zero argument) constructor
//        _hours = 0;
//        _mins = 0;
//        _secs = 0;}
//        
//    Time(int hours, int mins, int secs){ /* overloaded (3 arg) constructor */
//        _hours = hours;
//        _mins = mins;
//        _secs = secs;}
//        
//    ~Time() // default destructor
//        {}
//
//    void display(){
//        cout << _hours << ":" << _mins << ":" << _secs;
//        cout << endl;}
//
//    void add_time(Time t1, Time t2){
//        _secs = (t1._secs + t2._secs) % 60;
//        _mins = ((t1._mins + t2._mins) + (t1._secs + t2._secs) / 60) % 60;
//        _hours = ((t1._hours + t2._hours) +
//        (t1._mins + t2._mins + (t1._secs + t2._secs) / 60) / 60) % 24;}
//
//    friend Time operator+(Time t1, Time t2){
//        ushort secs = (t1._secs + t2._secs) % 60;
//        ushort mins = ((t1._mins + t2._mins) + (t1._secs + t2._secs) / 60) % 60;
//        ushort hours = ((t1._hours + t2._hours) +
//        (t1._mins + t2._mins + (t1._secs + t2._secs) / 60) / 60) % 24;
//        
//        Time temp(hours, mins, secs);
//        
//        return temp;}
//
//    friend ostream& operator<<(std::ostream& cout, Time& t){
//        cout << t._hours << ":" << t._mins << ":" << t._secs;
//        return cout;}
//};
//
//
//int main(){
//    Time time1(8, 45, 44);
//    Time time2(14, 5, 57);
//    
//    Time time3(0, 0, 0);
//    Time time4(0, 0, 0);
//    time1.display();
//    cout << endl;
//    
//    time3.add_time(time1, time2);/* time3 calls add function and is assigned time1+time2 */
//    time3.display(); //output
//    time4 = time1 + time2;
//    cout<<time4; /* using the overloaded operator + to directly add times
//    and then the overloaded console output statement to the screen*/
//    cout << endl;
//    return 0;
//}


// Declaring arrays of time objects
//****************************
//Time time4; /* note that this object upon creation will be set automatically to (0:0:0)*/
//Time array[4]={time4}; /* here we declare an array of type Time to hold 4 objects and we have initialized it to time3, i.e. (0:0:0)*/

//Within this same code we can define a pointer to a time object and then, have the pointer point to an array which we created above.
//Time* ptr; // declare a pointer of Time class
//ptr=&array[0]; /* initialize pointer by pointing to base address of array
//*/
//for (int i=0; i<4; ++i){
//cout<<ptr[i]<<endl;}
///* the loop above now prints out the contents of the array given by ptr */
//We can also create dynamic variables of type Time. Consider the following member function for
//printing an object that has a pointer as part of the parameter list.
//Time* Time:: display(Time* ptr) /* function to display object to screen */
//{
//return ptr;}

// with the following code in the main body
//Time* ptr=new Time; /* create a dynamic variable to which ptr points */
//Time time1(8,45,44); // instantiate an object
//*ptr=time1; // the dynamic variable is initialised
//ptr->display(); /* this is how a pointer is connected to a member function – in this case it prints the variable */
//delete ptr; // delete ptr_time;
//ptr=NULL; // return null address
//We note a new form of syntax -> used when the pointer acted on the function. This is how a pointer is
//connected to a member function, and is an example of syntactic sugar.
//In module 2 we used the fstream library to access the class ofstream and hence create an object
//out which allowed us to output results to an Excel file (or .dat and .txt formats). This idea also extends
//to print out objects in various files, in precisely the same way.
//The next few classes are mathematical examples. They give an important insight to the construction
//and running of classes.




//***********************************************************************
//******************** BS Option pricing model **************************
//***********************************************************************

//const double pi=4.0*atan(1.0); //define constant pi=3.142
//
//class BS_Pricing{
//private:
//    int _stock, _strike;
//    float _vol, _int_rate, _divvy, _tau;
//
//public:
//    BS_Pricing(){ // default (zero argument) constructor
//        _stock = 0;
//        _strike = 0;
//        _vol = 0;
//        _int_rate = 0;
//        _divvy = 0;
//        _tau = 0;}
//
//    BS_Pricing(int stock, int strike, float vol, float int_rate,
//               float divvy, float tau){
//                _stock = stock;
//                _strike = strike;
//                _vol = vol;
//                _int_rate = int_rate;
//                _divvy = divvy;
//                _tau = tau;}
//        
//    ~BS_Pricing() // destructor
//        {}
//    
//    float CDF(double X){ // approximation
//        const double a1=0.319381530, a2=-0.356563782, a3=1.781477937, a4=-1.821255978, a5=1.330274429; // define all constants a(i)
//        double x=0, N, CDF, n, k=0;
//
//        x=abs(X);
//        k=1/(1+0.2316419*x);
//        n=(1/sqrt(2*pi))*exp(-0.5*x*x);
//        N=1-n*(a1*k+a2*k*k+a3*pow(k,3)+a4*pow(k,4)+a5*pow(k,5));
//        CDF=N;
//        if (X<0){ /* we calculate for +ve X and then use the symmetry property of the distribution to obtain the CDF for –ve values*/
//            CDF=1-N;
//        }
//        return CDF;
//    }
//    
//    float Vanilla_Option_value(int stock, int strike, float vol, float int_rate,
//                               float divvy, float tau, int is_call=1){
//        
//        float d1=0, d2=0;
//        
//        d1 = (log(stock/strike) + (int_rate + divvy + 0.5*pow(vol,2)) * tau)/(vol * sqrt(tau));
//        d2 = d1 - vol*sqrt(tau);
//        
//        if (is_call==1){
//            return stock*CDF(d1) - strike*exp(-int_rate*tau)*CDF(d2);
//        }
//        else
//            return strike*exp(-int_rate*tau)*CDF(-d2) - stock*CDF(-d1);
//        
//        return 0;
//    }
//    
//    
//    float Binary_Option_value(int stock, int strike, float vol, float int_rate, float divvy, float tau, int is_call=1){
//        
//        float d1=0, d2=0;
//        
//        d1 = (log(stock/strike) + (int_rate + divvy + 0.5*pow(vol,2)) * tau)/(vol * sqrt(tau));
//        d2 = d1 - vol*sqrt(tau);
//        
//        if (is_call==1){
//            if (stock>strike){
//                return 1;
//            }
//            else
//                return 0;
//        }
//        else
//            if (stock>strike){
//                return 0;
//            }
//            else
//                return 1;
//        
//        return 0;
//    }
//
//};
//

//int main(){
//    float call,put, b_call, b_put;
//    BS_Pricing option_value;
//    call = option_value.Vanilla_Option_value(100, 95, 0.2, 0.03, 0.01, 0.25, 1);
//    put = option_value.Vanilla_Option_value(100, 95, 0.2, 0.03, 0.01, 0.25, 0);
//    b_call = option_value.Binary_Option_value(100, 95, 0.2, 0.03, 0.01, 0.25, 1);
//    b_put = option_value.Binary_Option_value(100, 95, 0.2, 0.03, 0.01, 0.25, 0);
//    cout<<"Call: "<< call << " Put: "<< put <<endl;
//    cout<<"Binary Call: "<< b_call << " Put: "<< b_put <<endl;
//    return 0;
//}


// GBM and Vasicek simulations. Perform multiple sample paths and output to files such as csv
//*******************************************************************************

//#include <ctime>
//#include <fstream>
//#include <iostream>
//#include <iomanip>
//#include <cstring>
//#include <cmath>
//
//using namespace std;
//
//const float pi = 3.14152;
//const double OneOverRootTwoPi = 0.398942280401433;
//
//double NormalDensity(double x);
//double normal(double, double);
//
//
//class SDE_generator{
//private:
//    float _stock;
//    float _vol;
//    int _N;
//    float _IR;
//    float _kappa;
//    float _theta;
//    
//    
//public:
//    // Constructor
//    SDE_generator(){
//        _stock = 0;
//        _vol = 0;
//        _N = 0;
//        _IR = 0;
//        _kappa = 0;
//        _theta = 0;
//    }
//    
//    // Definition
//    SDE_generator(float stock, float vol, int N, float IR,  float kappa, float theta){
//        _stock = stock;
//        _vol = vol;
//        _N = N;
//        _IR = IR;
//        _kappa = kappa;
//        _theta = theta;
//    }
//    
//    // Destructor
//    ~SDE_generator(){}
//    
//    vector<float> GBM_Generator(float stock, float vol, int N, float IR){
//        
//        srand((unsigned)time(NULL));
//        vector<float> array(N, 0);
//        
//        float dt = 1.0 / N; // step-size for time
//        array[0] = stock;
//        
//        for (unsigned short int i = 1; i <= N; i++) {
//            double time = i * dt;
//            double dX = normal(1.0, 0.0) * sqrt(dt);
//            double dS = stock * (IR * dt + vol * dX);
//            stock += dS;
//            array[i] = stock;
//        }
//        
//        return array;
//    }
//    
//    vector<float> Vasicek_Generator(float vol, int N, float IR, float kappa, float theta){
//        
//        float current_rate = 1.0;
//        srand((unsigned)time(NULL));
//        
//        vector<float> array(N, 0);
//
//        float dt = 1.0 / N; // step-size for time
//        array[0] = current_rate;
//        
//        for (unsigned short int i = 1; i <= N; i++) {
////            double time = i * dt;
//            double dX = normal(1.0, 0.0) * sqrt(dt);
//            double dr = kappa*(theta - IR) * dt + vol * dX; // - Milstein approx
//            current_rate += dr;
//            array[i] = current_rate;
//        }
//        
//        return array;
//    }
//};
//
//int main(){
//    int stock, N;
//    float vol, kappa, theta, IR;
//    vector<float> gbm_output, vasi_output;
//    
//    cout.precision(5);
//    cout << "What is the parameters: " <<endl;
//    cin >> stock >> vol >> kappa >> theta >> N >> IR;
//    
//    SDE_generator sde_gen;
//    gbm_output = sde_gen.GBM_Generator(stock, vol, N, IR);
//    vasi_output = sde_gen.Vasicek_Generator(vol, N, IR, kappa, theta);
//    
//    ofstream out;
//    out.precision(5);
//    out.open("C++/Code/Code/GBM.xls");
//    
//    for (int i = 1; i <= N; i++){
//        out << gbm_output[i] << endl;
//    };
//    
//    out.close();
//    
//    ofstream vasi_out;
//    vasi_out.precision(5);
//    vasi_out.open("/Users/a44791/Documents/Uni/AppCompFin/C++/Code/Code/Vasicek.xls");
//    
//    for (int i = 1; i <= N; i++){
//        vasi_out << vasi_output[i] << endl;
//    };
//    
//    vasi_out.close();
//    
//    return 0;
//}
//    // while: Generate random numebrs, calc GBM/vasicek at each t, store in array and output to file
//
//
//double NormalDensity(double x){
//    return OneOverRootTwoPi*exp(-x*x/2);
//    }
//
//
////Polar Marsaglia Transformation
//double normal(double std, double mean){
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



