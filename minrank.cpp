#include "minrank.h"

auto initialize_minor_table(){
    std::array<size_t , N*N*N*N> table;

    size_t ctr = 0;
    for (size_t a = 0; a < N; a++)
    {
        for (size_t b = a+1; b < N; b++)
        {
            for (size_t c = b+1; c < N; c++)
            {
                for (size_t d = c+1; d < N; d++)
                {
                    table[a*N*N*N + b*N*N + c*N + d] = ctr++;
                }
            }
        }
    }
    return table;
}

size_t get_minor_label(int i0, int i1, int i2, int i3){
    const static auto table = initialize_minor_table();
    return table[i0*N*N*N + i1*N*N + i2*N + i3];
}

Matrix vec_to_Mat(size_t width, const Vector &vec){
    Matrix M(vec.size()/width, width);
    for (size_t i = 0; i < M.rows ; i++)
    {
        for (size_t j = 0; j < M.cols; j++)
        {
            MData(M,i,j) = vec[(i*M.cols)+j];
        }
    }
    return M;
}

bool is_zero(const Vector &v){
    for (const long long &vi : v)
    {
        if((vi%Q) != 0){
            return false;
        }
    }
    return true;
}

void remove_trailing_zeros(Vector &v){
    while( v.size() > 0 && ((v[v.size()-1]%Q)+Q)%Q == 0){
        v.resize(v.size()-1);
    }
}

// a = quot*b + rem
void poly_division(const Vector &a, const Vector &B, Vector &quot, Vector &rem){

    Vector b = B;

    remove_trailing_zeros(b);

    if(b[b.size()-1] != 1){
        //std::cout << "Leading terms coefficient != 1.";
        long long inv = inverse_mod_Q(b[b.size()-1]);
        for (size_t i = 0; i < b.size(); i++)
        {
            b[i] = (((b[i]*inv)%Q)+Q)%Q;
        }
    }

    rem = a;

    if(a.size() < b.size()){
        quot = Vector(0);
        return;
    }

    quot.resize(a.size()-b.size()+1);

    for (int i = ((int)a.size())-((int)b.size()); i >= 0; i--)
    {
        quot[i] = ((rem[i+b.size()-1]%Q) +Q )% Q;
        for (size_t j = 0; j < b.size(); j++)
        {
            rem[i+j] -= quot[i]*b[j];
            rem[i+j] = ((rem[i+j]%Q)+Q)%Q;
        }
    }

    remove_trailing_zeros(rem);
}

Vector poly_mul(const Vector &a, const Vector &b){
    Vector mul(a.size()+b.size()-1);
    for (size_t i = 0; i < a.size(); i++)
    {
        for (size_t j = 0; j < b.size(); j++)
        {
            mul[i+j] += a[i]*b[j];
        }
    }
    for (size_t i = 0; i < mul.size(); i++)
    {
        mul[i] %= Q;
    }
    return mul;
}

Vector poly_gcd(const Vector &a, const Vector &b){
    Vector X = a;
    Vector Y = b;

    while (true)
    {
        if(X.size() < Y.size()){
            std::swap(X,Y);
        }

        remove_trailing_zeros(Y);
        if(Y.size() == 0){
            return X;
        }

        Vector Quo;
        Vector newX;

        poly_division(X,Y, Quo, newX);
        X = newX;
    }
}

Vector poly_lcm(const Vector &a, const Vector &b){
    Vector mul = poly_mul(a,b);
    Vector gcd = poly_gcd(a,b);
    Vector rem,lcm;
    poly_division(mul,gcd,lcm,rem);
    assert(rem.size() == 0);
    return lcm;
}

Vector exp_mod_poly(const Vector &base, size_t exp, const Vector &modulus){
    Vector out(1,1);
    Vector X,Y;

    Vector temp = base;
    while(exp > 0){
        if(exp%2){
            out = poly_mul(out, temp);
            poly_division(out,modulus, X,Y);
            out = Y;
        }

        temp = poly_mul(temp,temp);
        poly_division(temp,modulus, X,Y);
        temp = Y;
        exp /= 2;
    }
    return out;
}

Vector min_poly(const Matrix M){
    assert(M.cols == M.rows);

    Vector minpoly(1,1);

    #define TRIES 5
    for (size_t tr = 0; tr < TRIES; tr++)
    {
        Vector L(M.cols), R(M.cols);
        randomize_vector(L);
        randomize_vector(R);

        Vector Seq(2*M.cols+1);
        for (size_t i = 0; i < 2*M.cols+1; i++)
        {
            Seq[i] = L*R;
            R = M*R;
        }
        
        Matrix RetardedBM(M.cols+1, M.cols+1);
        for (size_t i = 0; i < M.cols+1; i++)
        {
            for (size_t j = 0; j < M.cols+1; j++)
            {
                MData(RetardedBM,i,j) = Seq[i+j];
            }
        }
        
        RetardedBM.EF();
        size_t r = RetardedBM.rank_EF();

        Vector poly = RetardedBM.kernel_vector_EF(M.cols - r); 
        minpoly = poly_lcm(minpoly, poly);
    }       

    // compute gcd with X^Q-X
    Vector field_eq = {0,1};
    field_eq = exp_mod_poly(field_eq,Q, minpoly);
    field_eq[1] = (Q-1 + field_eq[1])%Q;
    return poly_gcd(field_eq, minpoly);
}

long long horner(Vector poly, long long x){
    long long eval = 0;
    for (int i = poly.size()-1 ; i >= 0; i--)
    {
        eval = (eval*x + poly[i])%Q;
    }
    return eval;
}

Vector find_lambdas(const Matrix &A, const Matrix &B){
    std::vector<long long >lambdas(1,Q);

    Matrix AA(A.cols,A.cols);
    Matrix BB(B.cols,B.cols);

    for (size_t i = 0; i < A.rows; i++)
    {
        long long r = random()%Q;
        AA.set_row(i%A.cols, AA.row(i%A.cols) + (r* A.row(i)));
        BB.set_row(i%A.cols, BB.row(i%A.cols) + (r* B.row(i)));
    }
    
    Matrix B_invA(A.cols,A.cols);
    int res = multiply_by_inverse(BB,AA,B_invA);

    if(res < 0){
        return Vector(0);
    }

    Matrix minusB_invA(A.cols,A.cols);
    minusB_invA -= B_invA;

    return min_poly(minusB_invA);
}

int salvage_attempt(Matrix &CM, const std::vector<Matrix> Matrices, std::vector<Vector> &sols){
    auto CMKernel = CM.kernel_EF();
    assert(CMKernel.size()<=6);

    Vector ea(Matrices.size());
    Vector eb(Matrices.size());
    randomize_vector(ea);
    randomize_vector(eb);

    // construct A and B
    Matrix A(N_CHOOSE_R,CMKernel.size());
    Matrix B(N_CHOOSE_R,CMKernel.size());
    for (size_t kvi=0; kvi < CMKernel.size(); kvi++)
    {
        Matrix KM = vec_to_Mat(Matrices.size(), CMKernel[kvi]);
        A.set_col(kvi, KM*ea);
        B.set_col(kvi, KM*eb);
    }

    #define SMALL_HEIGHT 16
    Matrix A_small(SMALL_HEIGHT,CMKernel.size());
    Matrix B_small(SMALL_HEIGHT,CMKernel.size());
    for (size_t i=0; i < SMALL_HEIGHT; i++)
    {
        A_small.set_row(i, A.row(i));
        B_small.set_row(i, B.row(i));
    }


    Vector mp = find_lambdas(A,B);
    size_t zeros = 0;
    // solve for lambda s.t. A + lambda*B is singular
    // TODO: do this in a less retarded way
    for (size_t lambda = 0; lambda <= Q; lambda++)
    {
        long long eval = horner(mp,lambda);
        if(eval != 0){
            continue;
        }
        zeros += 1;
        if(zeros == N+1){
            return -2;
        }

        Matrix C_small = A_small;
        add_mul(lambda,B_small,C_small);
        if(lambda == Q){
            C_small = B_small;
        }

        if(C_small.rank() == CMKernel.size()){
            continue;
        }

        Matrix C = A;
        add_mul(lambda,B,C);

        if(lambda == Q){
            C = B;
        }

        C.EF();

        size_t Crank = C.rank_EF();

        if(Crank == CMKernel.size()){
            continue;
        }

        if(Crank < CMKernel.size()-1){
            return -2;
        }

        Vector X = C.kernel_vector_EF();

        Vector BMV(N_CHOOSE_R*Matrices.size(), 0); 
        for (size_t i = 0; i < X.size(); i++)
        {
            BMV += ( X[i] * CMKernel[i] );
        }
        Matrix BM = vec_to_Mat(Matrices.size(), BMV);

        for (size_t i = 0; i < BM.rows; i++)
        {
            Vector R = BM.row(i);
            if (!is_zero(R)){

                Matrix CheckM(N,N);

                for (size_t j = 0; j < R.size(); j++)
                {
                    add_mul(R[j], Matrices[j],CheckM);
                }

                if(CheckM.rank() <= 4){
                    sols.push_back(R);
                }

                break;
            }
        }
    }
    
    return 0;
}

// return 0 if success, -1 if no solution, -2 if undecided
int MinRank_attempt(const std::vector<Matrix> Matrices, std::vector<Vector> &solutions){
    size_t number_of_matrices = Matrices.size();

    if(number_of_matrices == 1){
        size_t r = Matrices[0].rank();
        if (r <= 4){
            solutions.emplace_back(Vector(1,1));
            return 0;
        }
        return -1;
    }

//std::cout << "CM is  " << N*N_CHOOSE_R_PLUS_ONE << " by " <<N_CHOOSE_R*number_of_matrices << std::endl;
//std::cout << "number of matrices: " << number_of_matrices << std::endl;

    size_t CMrows = N_CHOOSE_R*number_of_matrices + 31;
    //size_t CMrows = N*N_CHOOSE_R_PLUS_ONE;
    Matrix CM(CMrows,N_CHOOSE_R*number_of_matrices);

    size_t row_ctr = 0;
    for (size_t i0 = 0; i0 < N; i0++)
    {
        for (size_t i1 = i0+1; i1 < N; i1++)
        {
            for (size_t i2 = i1+1; i2 < N; i2++)
            {
                for (size_t i3 = i2+1; i3 < N; i3++)
                {
                    for (size_t i4 = i3+1; i4 < N; i4++)
                    {
                        for (size_t row = 0; row < N; row++)
                        {
                            long long r = random()%Q;
                            for (size_t mat = 0; mat < number_of_matrices; mat++)
                            {
                                MData(CM,row_ctr,(get_minor_label(i1,i2,i3,i4)*number_of_matrices) + mat) += (r*MData(Matrices[mat],row,i0))%Q;
                                MData(CM,row_ctr,(get_minor_label(i0,i2,i3,i4)*number_of_matrices) + mat) += (r*(Q-MData(Matrices[mat],row,i1))%Q)%Q;
                                MData(CM,row_ctr,(get_minor_label(i0,i1,i3,i4)*number_of_matrices) + mat) += (r*MData(Matrices[mat],row,i2))%Q;
                                MData(CM,row_ctr,(get_minor_label(i0,i1,i2,i4)*number_of_matrices) + mat) += (r*(Q-MData(Matrices[mat],row,i3))%Q)%Q;
                                MData(CM,row_ctr,(get_minor_label(i0,i1,i2,i3)*number_of_matrices) + mat) += (r*MData(Matrices[mat],row,i4))%Q;
                            }
                            row_ctr = (row_ctr+1)% CMrows;
                        }
                    }
                }
            }
        }
    }
//TIC
    CM.EF();
//TOC(Echelon form)
    size_t r = CM.rank_EF();
    
//    std::cout << N_CHOOSE_R*number_of_matrices - r << std::endl; 

    if(r == N_CHOOSE_R*number_of_matrices){
        return -1;
    }

    if(r < N_CHOOSE_R*number_of_matrices-1){

        if(r >= N_CHOOSE_R*number_of_matrices-6){
            // try up to 3 times in case salvage is undecided
//TIC
            int ret = salvage_attempt(CM, Matrices, solutions);
//TOC(Salvage)
            if(ret == -2){
                solutions.resize(0);
                ret = salvage_attempt(CM, Matrices, solutions);

                if(ret == -2){
                std::cout << "retry third time: ";
                solutions.resize(0);
                ret = salvage_attempt(CM, Matrices, solutions);
                    if(ret == -2){
                        std::cout << " Failed.\n";
                    } else {
                        std::cout << " Ok.\n";
                    }
                }
            }
            return ret;
        }

        std::cout << "kernel size:" << N_CHOOSE_R*number_of_matrices - r << std::endl;
        return -2;
    }

    Vector kernel_vector(CM.kernel_vector_EF());

    for (size_t i = 0; i < N_CHOOSE_R; i++)
    {
        int non_zero = 0;
        for (size_t j = 0; j < number_of_matrices; j++)
        {
            if(kernel_vector[i*number_of_matrices + j] != 0){
                non_zero = 1;
                break;
            }
        }

        if(non_zero){
            Vector solution(number_of_matrices);
            for (size_t j = 0; j < number_of_matrices; j++)
            {
                solution[j] = kernel_vector[i*number_of_matrices + j];
            }

            Matrix SOLM(N,N);
            for (size_t i = 0; i < number_of_matrices; i++)
            {
                add_mul(solution[i],Matrices[i],SOLM);
            }

            if(SOLM.rank() <= 4){
                solutions.emplace_back(Vector(solution));
                return 0;
            }

            break;
        }
    }

    return -1;
}

std::vector<Vector> enumerate_low_rank(const std::vector<Matrix> input_matrices){
    std::vector<Vector> rank_4_nodes;
    Vector guesses(input_matrices.size());
    assert(input_matrices.size()>=3);
    size_t pos = input_matrices.size()-1;

    while (guesses[0] <= 1)
    {
        if (pos < 0){
            break;
        }

        guesses[pos] += 1;

        if(guesses[pos] == Q){
            pos -= 1;
            continue;
        }

        int all_zero = 1;
        for (size_t i = 0; i < pos; i++)
        {
            if(guesses[i] != 0){
                all_zero = 0;
                break;
            }
        }

        if(all_zero && guesses[pos]>1){
            pos -= 1;
            continue;
        }        

        //for (size_t i = 0; i <= pos; i++)
        //{
        //    std::cout << guesses[i] << " " ;
        //}
        //std::cout << std::endl;
        
        // prepare MinRank problem
        std::vector<Matrix> Matrices(input_matrices.size()-pos,Matrix(N,N));
        std::vector<Vector> solutions(0);

        for (size_t i = 0; i <= pos; i++)
        {
            add_mul(guesses[i],input_matrices[i],Matrices[0]);
        }
        Matrices[0].reduce();

        for (size_t i = pos+1; i < input_matrices.size(); i++)
        {
            Matrices[i-pos] = input_matrices[i];
        }

        // solve minrank problem
        int result = MinRank_attempt(Matrices,solutions);

        if(result == -2){
            pos += 1;
            guesses[pos] = -1;
            continue;
        }

        if(result == -1){
            continue;
        }

        for(Vector &solution : solutions){

            normalize_vector(solution);

            if(solution[0] == 0){
                continue;
            }

            rank_4_nodes.emplace_back(input_matrices.size());
            Vector &SOL = rank_4_nodes.back();
            for (size_t i = 0; i <= pos; i++)
            {
                SOL[i] = guesses[i];
            }

            for (size_t i = pos+1; i < input_matrices.size(); i++)
            {
                SOL[i] = solution[i-pos];
            }
        }
    }

    return rank_4_nodes;
}
