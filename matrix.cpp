#include "matrix.h"

long long inverse_mod_Q(long long a){
    //EGCD
    a = ((a%Q)+Q) %Q;
    long long b = Q;
    long long c = 1;
    long long d = 0;

    while (a > 1)
    {
        int64_t q = b/a;
        b = b-q*a;
        d = d-q*c;

        int64_t temp;
        temp = a;
        a = b;
        b = temp;

        temp = c;
        c = d;
        d = temp;
    }
    return (Q+(c%Q)) % Q;
}

Matrix::Matrix(size_t rows, size_t cols)
{
    this->rows = rows;
    this->cols = cols;
    col_length = ((cols+3)/4)*4;
    data = Vector(col_length*rows,0);
}

void Matrix::reduce(){
    for (size_t i = 0; i < rows*col_length ; i++)
    {
        data[i] = (data[i]%Q+Q)%Q;
    }
}

void Matrix::randomize(){
    for (size_t i = 0; i < rows*col_length ; i++)
    {
        data[i] = random()%Q;
    }
}

Matrix& Matrix::operator+=(const Matrix& rhs){
    for (size_t i = 0; i < rows*col_length ; i++)
    {
        data[i] += rhs.data[i];
    }
    return *this;
}

Matrix& Matrix::operator-=(const Matrix& rhs){
    for (size_t i = 0; i < rows*col_length ; i++)
    {
        data[i] -= rhs.data[i];
    }
    return *this;
}

void add_mul(long long scalar, const Matrix& summand, Matrix& accumulator){
    scalar %= Q;
    assert(summand.cols == accumulator.cols);
    assert(summand.rows == accumulator.rows);
    for (size_t i = 0; i < summand.rows*summand.col_length ; i++)
    {
        accumulator.data[i] += scalar*summand.data[i];
    }
}

void Matrix::zero(){
    for (size_t i = 0; i < rows*col_length ; i++)
    {
        data[i] = 0;
    }
}

void Matrix::identity(){
    zero();
    for (size_t i = 0; i < std::min(rows,cols); i++)
    {
        Data(i,i) = 1;
    }
}

std::ostream& operator<<(std::ostream& os, const Matrix &Mat)
{
    for (size_t i = 0; i < Mat.rows; i++)
    {
        for (size_t j = 0; j < Mat.cols; j++)
        {
            os << MData(Mat,i,j) << " ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Vector &vec)
{
    for (long long vi : vec)
    {
        os << vi << " ";
    }
    return os;
}

void Matrix::EF()
{
    size_t col = 0;
    size_t row = 0;
    
    while(col < cols){
        // find nonzero entry in col
        size_t i = row;
        while (i < rows && Data(i,col)%Q == 0)
        {
            i ++;
        }

        if(i==rows){
            // no nonzero entry, move to next col
            col ++;
            continue;
        }

        if(i > row){
            // swap
            for (size_t j = col; j < cols; j++)
            {
                std::swap(Data(row, j),Data(i,j));
            }
        }

        long long inv = inverse_mod_Q(Data(row,col));
        for(size_t j=col; j< cols ; j++){
            //Data(row,j) = ((((Data(row,j)%Q)*inv) % Q) + Q) % Q;
            Data(row,j) = ((Data(row,j)%Q)*inv) % Q;
        }
        for (size_t i = row+1; i < rows; i++)
        {
            Data(i,col) %= Q;
            if(Data(i,col) != 0){
                long long scalar = (Q - (Data(i,col)%Q))%Q;
                for (size_t j = col+1; j < cols; j++)
                {
                    Data(i,j)   += Data(row,j)*scalar;
                }
            }
            Data(i,col) = 0;
        }
        col ++;
        row ++;
    }
}

void Matrix::RREF()
{
    size_t col = 0;
    size_t row = 0;
    
    while(col < cols){
        // find nonzero entry in col
        size_t i = row;
        while (i < rows && Data(i,col)%Q == 0)
        {
            i ++;
        }

        if(i==rows){
            // no nonzero entry, move to next col
            col ++;
            continue;
        }

        if(i > row){
            // swap
            for (size_t j = col; j < cols; j++)
            {
                std::swap(Data(row, j),Data(i,j));
            }
        }

        long long inv = inverse_mod_Q(Data(row,col));

        for(size_t j=col; j< cols ; j++){
            Data(row,j) = ((((Data(row,j)%Q)*inv) % Q) + Q) % Q;
        }

        for (size_t i = 0; i < rows; i++)
        {
            if(i==row){
                continue;
            }

            Data(i,col) %= Q;
            if(Data(i,col) != 0){
                for (size_t j = col+1; j < cols; j++)
                {
                    Data(i,j) -= Data(row,j)*Data(i,col);
                }
            }
            Data(i,col) = 0;
        }

        col ++;
        row ++;
    }
}

size_t Matrix::rank_EF() const{
    size_t i = 0;
    size_t j = 0;
    while (j<rows){
        if(((Data(i,j)%Q)+Q)%Q == 1){
            i += 1;
        }
        j+= 1;
    }
    return i;
}

size_t Matrix::rank() const{
    Matrix MM(*this);
    MM.EF();
    return MM.rank_EF();
}

void randomize_vector(Vector &vec){
    for (long long &vi : vec)
    {
        vi = rand() % Q;
    }
}

template<size_t LEN> 
void print_array(const std::array<long long, LEN> array){
    for (size_t i = 0; i < LEN; i++)
    {
        std::cout << (array[i]%Q + Q)%Q << ' ';
    }
    std::cout << std::endl;
}
template void print_array<FORM_LEN>(const std::array<long long, FORM_LEN> array);

template<size_t LEN> 
int normalize_array(std::array<long long, LEN> &array){
    size_t i=0;
    while (i < LEN && array[i]%Q == 0)
    {
        array[i] = 0;
        i++;
    }

    if(i == LEN){
        return -1;
    }

    long long inv = inverse_mod_Q(array[i]);

    array[i] = 1;
    i++;
    while (i < LEN)
    {
        array[i] = (array[i]*inv)%Q;
        i++;
    }
    return 0;
}
template int normalize_array<FORM_LEN>(std::array<long long, FORM_LEN> &array);


int normalize_vector(Vector &v){
    size_t i=0;
    while (i < v.size() && v[i]%Q == 0)
    {
        v[i] = 0;
        i++;
    }

    if(i == v.size()){
        return -1;
    }

    long long inv = inverse_mod_Q(v[i]);

    v[i] = 1;
    i++;
    while (i < v.size())
    {
        v[i] = (v[i]*inv)%Q;
        i++;
    }
    return 0;
}

// returns -1 if there is no pivot in this collumn
// returns the row index otherwise
const int Matrix::pivot_row(const size_t col){
    int i = rows-1;
    while(i >= 0 && Data(i,col)%Q == 0)
    {
        i--;
    }
    if(i == -1 || ((Data(i,col)%Q)+Q)%Q != 1){
        return -1;
    }

    for (size_t j = 0; j < col; j++)
    {
        if(Data(i,j)%Q != 0){
            return -1;
        }
    }
    
    return i;
}

const Vector Matrix::random_kernel_vector_EF(){
    Vector RHS(rows);
    Vector kernel_vec(cols);
    for (int col = cols-1; col >= 0; col--)
    {
        int r = pivot_row(col);
        if(r == -1){
            kernel_vec[col] = random()%Q;
        }
        else{
            kernel_vec[col] = RHS[r];
        }

        for (size_t i = 0; i < rows; i++)
        {
            RHS[i] = (((RHS[i] - kernel_vec[col]*Data(i,col))%Q)+Q)%Q;
        }
    }
    return kernel_vec;
}

const Vector Matrix::kernel_vector_EF(size_t i){
    Vector RHS(rows);
    Vector kernel_vec(cols);
    for (int col = cols-1; col >= 0; col--)
    {
        int r = pivot_row(col);
        if(r == -1){
            if(i == 0){
                kernel_vec[col] = 1;
            }
            else{
                kernel_vec[col] = 0;
            }
            i--;
        }
        else{
            kernel_vec[col] = RHS[r];
        }

        for (size_t i = 0; i < rows; i++)
        {
            RHS[i] = (((RHS[i] - kernel_vec[col]*Data(i,col))%Q)+Q)%Q;
        }
    }
    return kernel_vec;
}

const std::vector<Vector> Matrix::kernel_EF(){
    size_t rank = rank_EF();

    std::vector<Vector> kernel(cols-rank);
    for (size_t i = 0; i < cols-rank; i++)
    {
        kernel[i] = kernel_vector_EF(i);
    }
    return kernel;
}

const std::vector<Vector> Matrix::kernel(){
    Matrix EF = *this;
    EF.EF();
    return EF.kernel_EF();
}

Vector Matrix::operator*(const Vector& vec) const{
    assert(vec.size() == cols);

    Vector out(rows);
    for (size_t i = 0; i < rows; i++)
    {
        out[i] = 0;
        for (size_t j = 0; j < cols; j++)
        {
            out[i] += Data(i,j)*vec[j];
        }
        out[i] %= Q;
        out[i] += Q;
        out[i] %= Q;
    }
    return out;    
}

int multiply_by_inverse(const Matrix &mat, const Matrix &in, Matrix &out){
    size_t size = mat.cols;
    assert(mat.rows == size);
    assert(in.rows == size);
    assert(in.cols == size);
    assert(out.rows == size);
    assert(out.cols == size);

    Matrix temp(size,2*size);
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            MData(temp,i,j) = MData(mat,i,j);
        }
    }

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            MData(temp,i,size+j) = MData(in,i,j);
        }
    }

    temp.RREF();

    for (size_t i = 0; i < size; i++)
    {
        if(((MData(temp,i,i)%Q)+Q)%Q != 1){
           return -1; 
        }
    }

    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < size; j++)
        {
            MData(out,i,j) = MData(temp,i,size+j);
        }
    }

    return 0;
}

long long operator*(const Vector &v1, const Vector &v2){
    assert(v1.size() == v2.size());
    long long ans = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
        ans += v1[i]*v2[i];
    }
    ans %= Q;
    ans += Q;
    ans %= Q;
    return ans;
}

Vector operator*(long long c, const Vector &v){
    Vector multiple(v);
    c = (c+Q)%Q;
    for (long long &mi : multiple)
    {
        mi = (mi*c) % Q;
    }
    return multiple;
}

Vector& operator+=(Vector &v_acc, const Vector &v){
    assert(v_acc.size() == v.size());
    for (size_t i = 0; i < v_acc.size(); i++)
    {
        v_acc[i] += v[i];
    }
    return v_acc;
}

Vector operator+(const Vector &v1, const Vector &v2){
    Vector Out = v1;
    Out += v2;
    return Out;
}

bool operator==(const Vector &v1, const Vector &v2){
    if(v1.size() != v2.size()){
        return false;
    }
    for (size_t i = 0; i < v1.size(); i++)
    {
        if((v1[i] - v2[i])%Q != 0){
            return false;
        }
    }
    return true;
}

Vector Matrix::row(size_t r) const{
    Vector row(cols);
    for (size_t i = 0; i < cols; i++)
    {
        row[i] = Data(r,i);
    }
    return row;
}

Vector Matrix::col(size_t c) const{
    Vector col(rows);
    for (size_t i = 0; i < rows; i++)
    {
        col[i] = Data(i,c);
    }
    return col;
}

void Matrix::set_row(size_t r, const Vector &v){
    assert(v.size() == cols);

    for (size_t i = 0; i < cols; i++)
    {
        Data(r,i) = v[i];
    }
}

void Matrix::set_col(size_t c, const Vector &v){
    assert(v.size() == rows);

    for (size_t i = 0; i < rows; i++)
    {
        Data(i,c) = v[i];
    }
}

Matrix random_invertible_matrix(size_t size){
    Matrix M(size,size);
    while (true)
    {
        M.randomize();
        if(M.rank() == size){
            break;    
        }
    }
    return M;
}