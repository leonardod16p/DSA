#include <iostream>
#include <stdlib.h>
#include <malloc.h>
#include <typeinfo>

size_t portable_ish_malloced_size(const void* p) {
    return malloc_usable_size((void*)p);
}

using namespace std;

/*
Classes Definitions
*/

class c_matrix {
    int m, n;
    double* matrix;

public:

    int getfirstdimension();
    int getseconddimension();
    int setfirstdimension(int a);


    double* getmatrix();
    void setmatrix(double value, double* address);

    c_matrix(int a, int b) {
        m = a;
        n = b;
        //calculate the memory allocation
        unsigned int byte_size_vec = m * n * sizeof(double);
        if (!(matrix = (double*)malloc(byte_size_vec))) {
            cout << "No memory" << endl;
            exit(1);
        }
        //cout << "First Element Address: " << matrix << endl;
        //cout << "Last Element Address: " << matrix + m * n - 1 << endl;
        //cout << "Sizeof: " << byte_size_vec << endl;
        size_t true_length = portable_ish_malloced_size(matrix);
        //printf("%zu\n", true_length);
    };
    void show_matrix();
    void define_values();
    void null_matrix();
    c_matrix transpose();

    c_matrix operator+(c_matrix const& obj)
    {
        c_matrix result(obj.m, obj.n);
        int line = 0;
        int column = 0;
        int i = 0;
        while (i < m) {
            column = 0;
            line = i * n;
            while (column < n) {
                *(result.matrix + line + column) = *(matrix + line + column) + *(obj.matrix + line + column);
                ++column;
            }
            ++i;
        }
        return result;
    }

    friend c_matrix operator*(double, c_matrix);
    friend c_matrix operator*(c_matrix, double);

    c_matrix operator*(c_matrix const& obj) {
        //Two matrices multiplication of size m*n and n*p
        //Resulting matrix must be of m*p size
        //m x n * n x p
        int p = obj.n;
        if (n == obj.m) {
            c_matrix result(m, p); //deveria ser m e nao obj.m
            result.null_matrix();
            int a_line = 0;
            int b_line = 0;
            int c_line = 0;
            double sum = 0;
            double a_ik = 0;
            double b_kj = 0;
            for (int i = 0; i < m; ++i) {
                a_line = i * n;
                for (int j = 0; j < p; ++j) {
                    sum = 0;
                    for (int k = 0; k < n; ++k) {
                        b_line = k * p;
                        a_ik = *(matrix + a_line + k);
                        b_kj = *(obj.matrix + b_line + j);
                        sum = sum + a_ik * b_kj;
                    }
                    *(result.matrix + i*p + j) = sum;
                }
            }
            return result;
        }
        else {
            cout << "Matrix multiplication Undefined" << endl;
            cout << "First Matrix" << endl;
            cout << "m:  " << m << "\t" << "n:  " << n << "\t" << endl;
            cout << "Second Matrix" << endl;
            cout << "m:  " << obj.m << "\t" << "n:  " << obj.n << "\t" << endl;
            throw "You cannot perform this operation!" <<;
        }
    }

};

class c_square_matrix : public c_matrix {
    int m;
    double* matrix;
    
public:

    int getfirstdimension();
    int getseconddimension();
    using c_matrix::c_matrix;
    friend c_square_matrix operator*(double, c_square_matrix);
    friend c_square_matrix operator*(c_square_matrix, double);
    c_square_matrix(int a) :c_matrix{a,a} {
        m = a;
        //calculate the memory allocation
        unsigned int byte_size_vec = m * m * sizeof(double);
        if (!(matrix = (double*)malloc(byte_size_vec))) {
            cout << "No memory" << endl;
            exit(1);
        }
        //cout << "First Element Address: " << matrix << endl;
        //cout << "Last Element Address: " << matrix + m * m - 1 << endl;
        //cout << "Sizeof: " << byte_size_vec << endl;
        size_t true_length = portable_ish_malloced_size(matrix);
        //printf("%zu\n", true_length);
    };
    
    double trace();

};

class c_upper_triangular_matrix : public c_square_matrix {
    int m;
    double* matrix;
public:

    int getfirstdimension();
    int getseconddimension();
    c_square_matrix converting();
    double* getmatrix();

    c_upper_triangular_matrix(int a) : c_square_matrix{ a,a }
    {
        m = a;
        //calculate the memory allocation
        //The same memory size as the lower
        int count = 0;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j <= i; ++j) {
                count = count + 1;
            }
        }
        unsigned int byte_size_vec = count * sizeof(double);
        if (!(matrix = (double*)malloc(byte_size_vec))) {
            cout << "No memory" << endl;
            exit(1);
        }
        //cout << "First Element Address: " << matrix << endl;
        //cout << "Last Element Address: " << matrix + m * m - 1 << endl;
        //cout << "Sizeof: " << byte_size_vec << endl;
        size_t true_length = portable_ish_malloced_size(matrix);
        //printf("%zu\n", true_length);
    };

    c_upper_triangular_matrix operator+(c_upper_triangular_matrix const& obj) {
        c_upper_triangular_matrix result(m);
        int line = 0;
        int index = 0;
        int discount = 0;
        for (int i = 1; i - 1 < m; i++) {
            for (int j = 0; j < discount + m; j++) {
                *(result.matrix + line + j) = *(matrix + line + j) + *(obj.matrix + line + j);
            }
            line = line + discount + m;
            discount--;
            cout << endl;
        }
        return result;
    }



    friend c_upper_triangular_matrix operator*(double, c_upper_triangular_matrix);
    friend c_upper_triangular_matrix operator*(c_upper_triangular_matrix, double);
    void define_values();
    void show_matrix();
    c_matrix transpose(c_upper_triangular_matrix);
    double det();
    double trace();
};

class c_lower_triangular_matrix : public c_square_matrix {
    int m;
    double* matrix;
public:
    int getfirstdimension();
    int getseconddimension();
    c_square_matrix converting();
    double* getmatrix();

    c_lower_triangular_matrix(int a) : c_square_matrix{ a } {
        m = a;
        //calculate the memory allocation
        int count = 0;
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j <= i; ++j) {
                count = count + 1;
            }
        }
        unsigned int byte_size_vec = count * sizeof(double);
        if (!(matrix = (double*)malloc(byte_size_vec))) {
            cout << "No memory" << endl;
            exit(1);
        }
        //cout << "First Element Address: " << matrix << endl;
        //cout << "Last Element Address: " << matrix + m * m - 1 << endl;
        //cout << "Sizeof: " << byte_size_vec << endl;
        size_t true_length = portable_ish_malloced_size(matrix);
        //printf("%zu\n", true_length);
    };

    c_lower_triangular_matrix operator+(c_lower_triangular_matrix const& obj) {
        c_lower_triangular_matrix result(m);
        int line = 0;
        int index = 0;

        for (int i = 1; i - 1 < m; i++) {
            line = (i - 1) * (i) / 2;
            for (int j = 0; j < i; j++) {
                index = line + j;
                *(result.matrix + index) = *(matrix + index) + *(obj.matrix + index);
            }
            cout << endl;
        }
        return result;
    }

    c_lower_triangular_matrix operator*(c_lower_triangular_matrix const& obj) {
        c_lower_triangular_matrix result(m);
        int a_line = 0;
        int b_line = 0;
        double sum = 0;
        double a_ik = 0;
        double b_kj = 0;
        int k = 0;

        for (int i = 1; i - 1 < m; ++i) { //iniciei i = 1 e deixei i-1
            a_line = (i - 1) * (i) / 2;
            for (int j = 0; j < i; ++j) {
                k = j;
                sum = 0;
                while  (k < i) {
                    b_line = (k) * (k+1) / 2;
                    a_ik = *(matrix + a_line + k);
                    b_kj = *(obj.matrix + b_line + j);
                    sum = sum + a_ik * b_kj;
                    ++k;
                }
                *(result.matrix + a_line + j) = sum;
            }
        }
        return result;
    }

    friend c_lower_triangular_matrix operator*(double, c_lower_triangular_matrix);
    friend c_lower_triangular_matrix operator*(c_lower_triangular_matrix, double);

    void define_values();
    void show_matrix();
    c_matrix transpose(c_lower_triangular_matrix);
    double det();
    double trace();
};

class c_diagonal_matrix : public c_square_matrix {
    int m;
    double* matrix;
public:
    c_square_matrix converting();
    c_diagonal_matrix(int a):c_square_matrix {a}
    {
        m = a;
        //calculate the memory allocation
        //The same memory size as the lower
        unsigned int byte_size_vec = m * sizeof(double);
        if (!(matrix = (double*)malloc(byte_size_vec))) {
            cout << "No memory" << endl;
            exit(1);
        }
        //cout << "Value of byte_size_vec: " << byte_size_vec << endl;
        //cout << "First Element Address: " << matrix << endl;
        //cout << "Last Element Address: " << matrix + m << endl;
        //cout << "Sizeof: " << byte_size_vec << endl;
        size_t true_length = portable_ish_malloced_size(matrix);
        //printf("%zu\n", true_length);
    };
    int getfirstdimension();
    int getseconddimension();
    double* getmatrix();
    void setmatrix(double value, double* address);
    
    friend c_diagonal_matrix operator*(double, c_diagonal_matrix);
    friend c_diagonal_matrix operator*(c_diagonal_matrix, double);
    void define_values();
    void show_matrix();
    double det();
    double trace();
    void null_matrix();

};

/*
c_matrix Methods
*/

int c_matrix :: getfirstdimension(){
    return m;
}

int c_matrix :: getseconddimension(){
    return m;
}

double* c_matrix :: getmatrix(){
    return matrix;
}

void c_matrix :: setmatrix(double value, double* address){
    *address = value;
};

void c_matrix::define_values() {
    cout << "To enter an element, type the number you want and press enter: " << endl;
    int line = 0;
    int column = 0;
    int i = 0;
    while (i < m) {
        column = 0;
        line = i * n;
        while (column < n) {
            cin >> *(matrix + line + column);
            ++column;
        }
        ++i;

    }
}

void c_matrix::show_matrix() {
    int line = 0;
    int column = 0;
    int i = 0;
    while (i < m) {
        column = 0;
        line = i * n;
        while (column < n) {
            cout << *(matrix + line + column) << "\t";
            ++column;
        }
        cout << endl;
        ++i;

    }
}

void c_matrix::null_matrix() {
    int line = 0;
    int column = 0;
    int i = 0;
    while (i < m) {
        column = 0;
        line = i * n;
        while (column < n) {
            *(matrix + line + column) = 0;
            ++column;
        }
        ++i;

    }
}

c_matrix c_matrix::transpose() {
    c_matrix T(n, m);
    T.null_matrix();
    int a_line = 0;
    int b_line = 0;
    int i = 0;
    int j = 0;
    while (i < m) {
        a_line = i * n;
        j = 0;
        while (j < n) {
            b_line = j * m;
            *(T.matrix + b_line + i) = *(matrix + a_line + j);
            ++j;
        }
        ++i;
    }
    return T;
}

/*
c_square_matrix Methods
*/

int c_square_matrix :: getfirstdimension(){
    return m;
}

int c_square_matrix :: getseconddimension(){
    return m;
}

double c_square_matrix :: trace() {
    int line = 0;
    int column = 0;
    double sum = 0;
    while (column < m) {
        line = column * m;
        sum = sum + *(matrix + line + column);
        ++column;
    }
    return sum;
}


/*
Upper Matrix Methods
*/


double* c_upper_triangular_matrix :: getmatrix(){
    return matrix;
}

int c_upper_triangular_matrix :: getfirstdimension(){
    return m;
}

int c_upper_triangular_matrix :: getseconddimension(){
    return m;
}

void c_upper_triangular_matrix::define_values() {
    cout << "To enter an element, type the number you want and press enter: " << endl;
    int line = 0;
    int index = 0;
    int discount = 0;
    for (int i = 1; i - 1 < m; i++) {
        for (int j = 0; j < discount + m; j++) {
            cin >> *(matrix + line + j);
        }
        line = line + discount + m;
        discount--;
        cout << endl;
    }

}

void c_upper_triangular_matrix::show_matrix() {
    //the formatting would be good only if all the elements of the matrix has size of 1 (in the sense of representation by strings, e.g size(142) = 3)
    //We could implement a better solution by creating an algorithm that checks the maximum size of the elements of the matrix
    //Then, we could use this maximum size to fix the amount of spaces needed
    int line = 0;
    int index = 0;
    int discount = 0;
    string space = "";
    for (int i = 1; i - 1 < m; i++) {
        cout << space;
        for (int j = 0; j < discount + m; j++) {
            cout <<  *(matrix + line + j) << " ";
        }
        space = space + " " + " ";
        line = line + discount + m;
        discount--;
        cout << endl;
    }
}

double c_upper_triangular_matrix::det() {
    int line = 0;
    int index = 0;
    int discount = 0;
    double prod = 1;
    for (int i = 1; i - 1 < m; i++) {
            prod = prod* (*(matrix + line ));
            line = line + discount + m;
            discount--;
        }
    cout << endl;
    return prod;
}

double c_upper_triangular_matrix::trace(){
    int line = 0;
    int index = 0;
    int discount = 0;
    double sum = 0;
    for (int i = 1; i - 1 < m; i++) {
            sum = sum* (*(matrix + line ));
            line = line + discount + m;
            discount--;
        }
    cout << endl;
    return sum;
}

c_matrix c_upper_triangular_matrix :: transpose(c_upper_triangular_matrix A) { 

    //Low Quality code
    c_matrix result(A.getfirstdimension(), A.getfirstdimension());

    result = A.converting();
    result = result.transpose();
    return result;
}

/*
Lower Matrix Methods
*/


double* c_lower_triangular_matrix :: getmatrix(){
    return matrix;
}


int c_lower_triangular_matrix :: getfirstdimension(){
    return m;
}

int c_lower_triangular_matrix :: getseconddimension(){
    return m;
}


void c_lower_triangular_matrix::define_values() {
    cout << "To enter an element, type the number you want and press enter: " << endl;
    int line = 0;
    int index = 0;

    for (int i = 1; i - 1 < m; i++) {
        line = (i - 1) * (i) / 2;
        for (int j = 0; j < i; j++) {
            index = line + j;
            cin >> *(matrix + index);
        }
        cout << endl;
    }
}

void c_lower_triangular_matrix::show_matrix(){
    int line = 0;
    int index = 0;

    for (int i = 1; i - 1 < m; i++) {
        line = (i - 1) * (i) / 2;
        for (int j = 0; j < i; j++) {
            index = line + j;
            cout <<  *(matrix + index) << " ";
        }
        cout << endl;
    }
}

double c_lower_triangular_matrix::det() {
    double prod = 1;
    int line = 0;
    int diagonal_index = 0;

    for (int i = 0; i < m; i++) {
        line = i * (i + 1) / 2;
        diagonal_index = line + i;
        prod = prod * (*(matrix + diagonal_index));
        }

    return prod;
}

double c_lower_triangular_matrix::trace(){
    double sum = 0;
    int line = 0;
    int diagonal_index = 0;

    for (int i = 0; i < m; i++) {
        line = i * (i + 1) / 2;
        diagonal_index = line + i;
        sum = sum + (*(matrix + diagonal_index));
        }

    return sum;
}


c_matrix c_lower_triangular_matrix::transpose(c_lower_triangular_matrix A) { 

    //Low Quality Code
    c_matrix result(A.getfirstdimension(), A.getfirstdimension());

    result = A.converting();
    result = result.transpose();
    return result;
 }

/*
Diagonal Matrix Methods
*/

int c_diagonal_matrix :: getfirstdimension(){
    return m;
}

int c_diagonal_matrix :: getseconddimension(){
    return m;
}

double* c_diagonal_matrix :: getmatrix(){
    return matrix;
}

void c_diagonal_matrix :: setmatrix(double value, double* address){
    *address = value;
};


void c_diagonal_matrix::define_values() {
    cout << "To enter an element, type the number you want and press enter: " << endl;
    int i = 0;
    while (i < m) {
        cin >> *(matrix + i);
        ++i;
    }
}

void c_diagonal_matrix::show_matrix() {
    int i = 0;
    while (i < m) {
        cout << *(matrix + i) << "\t";; 
        ++i;
    }
}

void c_diagonal_matrix::null_matrix() {
    int i = 0;
    while (i < m) {
        *(matrix + i) = 0;
        ++i;
    }
}

double c_diagonal_matrix::trace(){
    int i = 0;
    double sum = 0;
    while (i < m) {
        sum = sum + *(matrix + i);
        ++i;
    }
    return sum;
};

double c_diagonal_matrix::det(){
    int i = 0;
    double prod = 1;
    while (i < m) {
        prod = prod * *(matrix + i);
        ++i;
    }
    return prod;
}

/*
Methods to convert any type of matrix to c_square_matrix. We are going to use
this to overload operators on classes of different hierarchies
*/

c_square_matrix c_diagonal_matrix::converting(){
    double* matrix_result = nullptr;
    c_square_matrix result(m);
    result.null_matrix();
    matrix_result = result.getmatrix();
    int line = 0;
    int i = 0;
    while (i < m) {
        line = i * m;
        *(matrix_result + line + i) = *(matrix + i);
        ++i;
    }
    return result;
}
c_square_matrix c_lower_triangular_matrix::converting(){
    double* matrix_result = nullptr;
    c_square_matrix result(m);
    result.null_matrix();
    matrix_result = result.getmatrix();

    int a_line = 0;
    int b_line = 0;
    double value = 0;
    for (int i = 0; i  < m; i++) {
        a_line = i * m;
        b_line = (i + 1) * (i) / 2;
        for (int j = 0; j - 1 < i; j++) {
            *(matrix_result + a_line + j) = *(matrix + b_line + j);
        }
        cout << endl;
    }
    return result;
}
c_square_matrix c_upper_triangular_matrix::converting(){
    double* matrix_result = nullptr;
    c_square_matrix result(m);
    result.null_matrix();
    matrix_result = result.getmatrix();
    int a_line = 0;
    double value = 0;

    int i = 0;
    int j = 0;
    int k = 0;
    while(i < m) {
        a_line = (i) * m;
        while(j < m){
            *(matrix_result + a_line + j) = *(matrix + k);
            j++;
            k++;
        }
        i++;
        j = i;
        cout << endl;
    }

    return result;
}

/*
Defining Operations
*/

c_matrix operator*(double lambda, c_matrix A){
    int line = 0;
    int column = 0;
    int i = 0;
    int m = A.m;
    int n = A.n;
    c_matrix result(m,n);
    double* matrix = nullptr;
    matrix = A.matrix;
    while (i < m) {
        column = 0;
        line = i * n;
        while (column < n) {
            *(result.matrix + line + column) = lambda*(*(matrix + line + column)) ;
            ++column;
        }
        ++i;
    }
    return result;
};

c_matrix operator*(c_matrix A, double lambda){
    int line = 0;
    int column = 0;
    int i = 0;
    int m = A.m;
    int n = A.n;
    c_matrix result(m,n);
    double* matrix = nullptr;
    matrix = A.matrix;
    while (i < m) {
        column = 0;
        line = i * n;
        while (column < n) {
            *(result.matrix + line + column) = lambda*(*(matrix + line + column)) ;
            ++column;
        }
        ++i;
    }
    return result;
};

c_upper_triangular_matrix operator*(c_upper_triangular_matrix obj1, c_upper_triangular_matrix obj2 ) {
    int dimension = obj1.getfirstdimension();
    c_square_matrix A(dimension);
    A = obj1.converting();

    c_square_matrix B(dimension);
    B = obj2.converting();

    c_matrix C(dimension,dimension);
    C = (A*B);
    double* transform_address = nullptr; 
    double* final_matrix = nullptr; 
    c_upper_triangular_matrix result(dimension);
    transform_address = C.getmatrix();
    final_matrix = result.getmatrix();
    int a_line = 0;
    double value = 0;

    int i = 0;
    int j = 0;
    int k = 0;
    while(i < dimension) {
        a_line = (i) * dimension;
        while(j < dimension){
            result.setmatrix(*(transform_address + a_line + j),final_matrix + k);
            j++;
            k++;
        }
        i++;
        j = i;
        cout << endl;
    }


    return result;
    }

c_square_matrix operator*(double lambda, c_square_matrix A){
    int line = 0;
    int column = 0;
    int i = 0;
    int m = A.m;
    double value = 0;
    double* matrix = nullptr;
    double* result_matrix = nullptr;

    c_square_matrix result(m);


    result_matrix = result.getmatrix();
    matrix = A.getmatrix();
    while (i < m) {
        column = 0;
        line = i * m;
        while (column < m) {
            value = lambda*(*(matrix + line + column));
            result.setmatrix(value,result_matrix + line + column);
            ++column;
        }
        ++i;
    }
    return result;
};

c_square_matrix operator*(c_square_matrix A, double lambda){
    int line = 0;
    int column = 0;
    int i = 0;
    int m = A.m;
    double value = 0;
    double* matrix = nullptr;
    double* result_matrix = nullptr;

    c_square_matrix result(m);


    result_matrix = result.getmatrix();
    matrix = A.getmatrix();
    while (i < m) {
        column = 0;
        line = i * m;
        while (column < m) {
            value = lambda*(*(matrix + line + column));
            result.setmatrix(value,result_matrix + line + column);
            ++column;
        }
        ++i;
    }
    return result;
};

c_upper_triangular_matrix operator*(double lambda, c_upper_triangular_matrix A){
    int line = 0;
    int index = 0;
    int discount = 0;
    c_upper_triangular_matrix result(A.m);
    for (int i = 1; i - 1 < A.m; i++) {
        for (int j = 0; j < discount + A.m; j++) {
            *(result.matrix + line + j) = lambda*(*(A.matrix + line + j));
        }
        line = line + discount + A.m;
        discount--;
        cout << endl;
    }
    return result;
}

c_upper_triangular_matrix operator*(c_upper_triangular_matrix A, double lambda){
    int line = 0;
    int index = 0;
    int discount = 0;
    c_upper_triangular_matrix result(A.m);
    for (int i = 1; i - 1 < A.m; i++) {
        for (int j = 0; j < discount + A.m; j++) {
            *(result.matrix + line + j) = lambda*(*(A.matrix + line + j));
        }
        line = line + discount + A.m;
        discount--;
        cout << endl;
    }
    return result;
}

c_lower_triangular_matrix operator*(c_lower_triangular_matrix A, double lambda){
    int line = 0;
    int index = 0;

    c_lower_triangular_matrix result(A.m);
    for (int i = 1; i - 1 < A.m; i++) {
        line = (i - 1) * (i) / 2;
        for (int j = 0; j < i; j++) {
            index = line + j;
            *(result.matrix + index) = lambda*(*(A.matrix + index));
        }
        cout << endl;
    }
    return result;
}

c_lower_triangular_matrix operator*(double lambda, c_lower_triangular_matrix A){
    int line = 0;
    int index = 0;

    c_lower_triangular_matrix result(A.m);
    for (int i = 1; i - 1 < A.m; i++) {
        line = (i - 1) * (i) / 2;
        for (int j = 0; j < i; j++) {
            index = line + j;
            *(result.matrix + index) = lambda*(*(A.matrix + index));
        }
        cout << endl;
    }
    return result;
}


c_diagonal_matrix operator*(double lambda, c_diagonal_matrix A){
    int i = 0;
    c_diagonal_matrix result(A.m);
    while (i < A.m) {
        *(result.matrix + i) = lambda*(*(A.matrix + i));
        ++i;
    }
    return result;
}

c_diagonal_matrix operator*(c_diagonal_matrix A, double lambda){
    int i = 0;
    c_diagonal_matrix result(A.m);
    while (i < A.m) {
        *(result.matrix + i) = lambda*(*(A.matrix + i));
        ++i;
    }
    return result;
}


c_square_matrix operator+(c_square_matrix obj1, c_square_matrix obj2){
        double* matrix1 = nullptr;
        double* matrix2 = nullptr;
        double* matrix3 = nullptr;
        int m = obj1.getfirstdimension();
        c_square_matrix result(m);
        result.null_matrix();
        matrix1 = obj1.getmatrix();
        matrix2 = obj2.getmatrix();
        matrix3 = result.getmatrix();
        int a_line = 0;
        int b_line = 0;
        int c_line = 0;
        double sum = 0;
        double a_ik = 0;
        double b_kj = 0;
        for (int i = 0; i < m; ++i) {
            a_line = i * m;
            for (int j = 0; j < m; ++j)
                *(matrix3 + i * m + j) = *(matrix1 + a_line + j) + *(matrix2 + a_line + j);
        }
        return result;
    }
c_square_matrix operator*(c_square_matrix obj1, c_square_matrix obj2){
        double* matrix1 = nullptr;
        double* matrix2 = nullptr;
        double* matrix3 = nullptr;
        int m = obj1.getfirstdimension();
        c_square_matrix result(m);
        result.null_matrix();
        matrix1 = obj1.getmatrix();
        matrix2 = obj2.getmatrix();
        matrix3 = result.getmatrix();
        int a_line = 0;
        int b_line = 0;
        int c_line = 0;
        double sum = 0;
        double a_ik = 0;
        double b_kj = 0;
        for (int i = 0; i < m; ++i) {
            a_line = i * m;
            for (int j = 0; j < m; ++j) {
                sum = 0;
                for (int k = 0; k < m; ++k) {
                    b_line = k * m;
                    a_ik = *(matrix1 + a_line + k);
                    b_kj = *(matrix2 + b_line + j);
                    sum = sum + a_ik * b_kj;
                    }
                *(matrix3 + i * m + j) = sum;
            }

        }
        return result;
    }


c_square_matrix operator+(c_upper_triangular_matrix obj1, c_diagonal_matrix obj2){
    c_square_matrix A(obj2.getfirstdimension());
    A = obj1.converting();

    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();
    c_square_matrix C(obj2.getfirstdimension());

    /*
    Converting lower to square
    */
    C = A+B;
    return C;
}   
c_square_matrix operator+(c_square_matrix obj1, c_upper_triangular_matrix obj2){

    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    C = obj1+B;
    return C;

}
c_square_matrix operator+(c_square_matrix obj1, c_lower_triangular_matrix obj2) {
    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    C = obj1+B;
    return C;
}
c_square_matrix operator+(c_square_matrix obj1, c_diagonal_matrix obj2){
    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    C = obj1+B;
    return C;
    }
c_square_matrix operator+(c_upper_triangular_matrix obj1, c_lower_triangular_matrix obj2){
    c_square_matrix A(obj2.getfirstdimension());
    A = obj1.converting();

    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());

    /*
    Converting lower to square
    */
    C = A+B;
    return C;
}
c_lower_triangular_matrix operator+(c_lower_triangular_matrix obj1, c_diagonal_matrix obj2){
    int m = 0;

    double* matrix1 = nullptr;
    double* matrix2 = nullptr;
    double* matrix3 = nullptr;

    m = obj1.getfirstdimension();
    matrix1 = obj1.getmatrix();
    matrix2 = obj2.getmatrix();
    c_lower_triangular_matrix result(m);


    matrix3 = result.getmatrix();
    int a_line = 0;
    double value = 0;
    
    for (int i = 0; i  < (m + 1) * (m) / 2; i++) {
        obj1.setmatrix(*(matrix1 + i),matrix3 + i);
        }
        cout << endl;

    for (int i = 0; i  < m; i++) {
        a_line = (i + 1) * (i) / 2;
        value = *(matrix1 + a_line + i) + *(matrix2 + i);
        obj1.setmatrix(value,matrix3 + a_line + i);
        }
        cout << endl;
    return result;
    }
c_diagonal_matrix operator+(c_diagonal_matrix obj1, c_diagonal_matrix obj2){
    int m = 0;
    int n = 0;
    double* matrix1 = nullptr;
    double* matrix2 = nullptr;
    double* matrix3 = nullptr;
    m = obj1.getfirstdimension();
    matrix1 = obj1.getmatrix();
    matrix2 = obj2.getmatrix(); 
    c_diagonal_matrix result(m);

    matrix3 = result.getmatrix();
    double value = 0;
    int i = 0;
    while (i < m) {
        value = *(matrix1 + i) + *(matrix2 + i);
        obj1.setmatrix(value,matrix3 + i);
        ++i;
    }
    cout << endl;
    return result;
}



c_square_matrix operator*(c_square_matrix obj1, c_lower_triangular_matrix obj2) {

    c_square_matrix A(obj2.getfirstdimension(),obj2.getfirstdimension());
    A = obj2.converting();

    c_square_matrix B(obj2.getfirstdimension(),obj2.getfirstdimension());

    B = obj1*A;
    
    return B;




}
c_square_matrix operator*(c_square_matrix obj1, c_diagonal_matrix obj2){
    c_square_matrix A(obj2.getfirstdimension(),obj2.getfirstdimension());
    A = obj2.converting();

    c_square_matrix B(obj2.getfirstdimension(),obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    B = obj1*A;
    return B;
}
c_square_matrix operator*(c_lower_triangular_matrix obj1, c_diagonal_matrix obj2){
        c_square_matrix A(obj2.getfirstdimension(),obj2.getfirstdimension());
        A = obj1.converting();

        c_square_matrix B(obj2.getfirstdimension(),obj2.getfirstdimension());
        B = obj2.converting();
        c_square_matrix C(obj2.getfirstdimension(),obj2.getfirstdimension());
        /*
        Converting lower to square
        */
        C = A*B;
    return C;
}
c_square_matrix operator*(c_upper_triangular_matrix obj1, c_diagonal_matrix obj2){
    c_square_matrix A(obj2.getfirstdimension());
    A = obj1.converting();
    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    C = A*B;
    return C;
}
c_square_matrix operator*(c_square_matrix obj1, c_upper_triangular_matrix obj2){
    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();
    c_square_matrix C(obj2.getfirstdimension());

    /*
    Converting lower to square
    */
    C = obj1*B;
    return C;
}
c_square_matrix operator*(c_upper_triangular_matrix obj1, c_lower_triangular_matrix obj2){
    c_square_matrix A(obj2.getfirstdimension());
    A = obj1.converting();

    c_square_matrix B(obj2.getfirstdimension());
    B = obj2.converting();

    c_square_matrix C(obj2.getfirstdimension());
    /*
    Converting lower to square
    */
    C = A*B;
    return C;

}
c_diagonal_matrix operator*(c_diagonal_matrix obj1, c_diagonal_matrix obj2){
    int m = 0;
    int n = 0;
    double* matrix1 = nullptr;
    double* matrix2 = nullptr;
    double* matrix3 = nullptr;
    m = obj1.getfirstdimension();
    matrix1 = obj1.getmatrix();
    matrix2 = obj2.getmatrix(); 
    c_diagonal_matrix result(m);

    matrix3 = result.getmatrix();
    double value = 0;
    int i = 0;
    while (i < m) {
        value = *(matrix1 + i) * *(matrix2 + i);
        obj1.setmatrix(value,matrix3 + i);
        ++i;
    }
    cout << endl;
    return result;
}




/*-------------------------
      Global Variables
--------------------------*/

//c_matrix Matrices[100];     //Limite de 100 matrizes


int main()
{

    //Testing every overload operators 

    c_matrix A(2,3);
    c_matrix B(3,2);
    A.define_values();
    B.define_values();
    cout << "Testing c_matrix multiplication and multiplication by scalar: " << endl;
    A.show_matrix();
    B.show_matrix();
    (A*B).show_matrix();
    (3*A).show_matrix();
    (A*B+A*B).show_matrix();
    c_square_matrix C(3);
    C.define_values();
    cout << "Testing c_square_matrix multiplication and multiplication by scalar: " << endl;
    C.show_matrix();
    (C+C).show_matrix();
    (C*C).show_matrix();
    (3*C).show_matrix();

    c_upper_triangular_matrix E(3);
    E.define_values();
    cout << "Testing c_upper_triangular_matrix multiplication and multiplication by scalar: " << endl;
    E.show_matrix();
    (E+E).show_matrix();
    (E*E).show_matrix();
    (3*E).show_matrix();
    E.transpose(E).show_matrix();

    c_lower_triangular_matrix G(3);
    G.define_values();
    cout << "Testing c_lower_triangular_matrix multiplication and multiplication by scalar: " << endl;
    G.show_matrix();
    (G+G).show_matrix();
    (G*G).show_matrix();
    (3*G).show_matrix();

    c_diagonal_matrix I(3);
    I.define_values();
    cout << "Testing c_diagonal_matrix multiplication and multiplication by scalar: " << endl;
    I.show_matrix();

    (I+I).show_matrix();
    (I*I).show_matrix();
    (3*I).show_matrix();
}
