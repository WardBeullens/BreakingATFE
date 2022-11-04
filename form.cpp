#include "form.h"
#include <unordered_map>
#include <thread>

namespace std {
  template <>
  struct hash<Vector>
  {
    size_t operator()(const Vector& vec) const
    {
      size_t hash = vec.size();

      for(auto &i: vec){
        hash *= Q;
        hash += i;
      }
      return hash;
    }
  };
}

void randomize_form(Form &f){
    for (size_t i = 0; i < FORM_LEN; i++)
    {
        f[i] = random()%Q;
    }    
}

Matrix form_to_mat(const Form f, const std::vector<long long> &v){
    Matrix M(N,N);

    M.zero();
    size_t ctr = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i+1; j < N; j++)
        {
            for (size_t k = j+1; k < N; k++)
            {
                MData(M,i,j) += f[ctr]*v[k];
                MData(M,j,k) += f[ctr]*v[i];
                MData(M,k,i) += f[ctr]*v[j];
                MData(M,j,i) -= f[ctr]*v[k];
                MData(M,i,k) -= f[ctr]*v[j];
                MData(M,k,j) -= f[ctr]*v[i];
                ctr++;
            }
        }
    }
    M.reduce();
    return M;
}

void brute_force_rank(const Form f, const unsigned int r, std::vector<long long> &v){
    while(1){
        randomize_vector(v);
        Matrix M(form_to_mat(f,v));
        if(M.rank() == r){
            normalize_vector(v);
            return;
        }
    }
}


std::vector<std::vector<long long>> enumerate_low_rank(const Form &f){
    std::vector<Matrix> matrices(N);
    for (size_t i = 0; i < N; i++)
    {
        std::vector<long long> ei(N,0);
        ei[i] = 1;
        matrices[i] = form_to_mat(f,ei);
    }
    return enumerate_low_rank(matrices);
}


std::vector<std::vector<long long>> enumerate_low_rank_neighbours(const Form &f, const std::vector<long long> &v){
    Matrix M = form_to_mat(f,v);
    M.EF();
    size_t r = M.rank_EF();
    std::vector<Matrix> matrices(N-r);
    std::vector<std::vector<long long>> kernel_vectors(N-r);
    for (size_t i = 0; i < N-r; i++)
    {
        kernel_vectors[i] = M.kernel_vector_EF(i);
        matrices[i] = form_to_mat(f,kernel_vectors[i]);
    }

    size_t solutions = 0;
    auto lrn = enumerate_low_rank(matrices);
    std::vector<std::vector<long long>> output(lrn.size(),std::vector<long long>(N));
    for( auto &s : lrn){
        for (size_t i = 0; i < s.size() ; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                output[solutions][j] += s[i]*kernel_vectors[i][j];
                output[solutions][j] %= Q;
                output[solutions][j] += Q;
                output[solutions][j] %= Q;
            }
            
        }
        normalize_vector(output[solutions]);
        solutions ++;
    }

    return output;
}

Vector random_neighbour(const Form &f, const Vector &v){
    while(true){

        Matrix S = random_invertible_matrix(N);
        Form randomized_f = act_on_form(f, S);
        
        Matrix S_inv(N,N);
        S_inv.identity();
        multiply_by_inverse(S,S_inv,S_inv);
        
        Vector randomized_v = S_inv * v;

        Matrix M = form_to_mat(randomized_f,randomized_v);
        M.EF();
        size_t r = M.rank_EF();
        assert(r == 4);
        std::vector<Matrix> matrices(N-r-1);
        std::vector<Vector> kernel_vectors = M.kernel_EF();

        for (size_t i = 0; i < kernel_vectors.size()-1; i++)
        {
            matrices[i] = form_to_mat(randomized_f, kernel_vectors[i]);
        }
        
        auto r4n = enumerate_low_rank(matrices);
        for( auto &s : r4n){
            Vector output(N,0);

            for (size_t i = 0; i < s.size(); i++)
            {
                output += s[i] * kernel_vectors[i];
            }

            output = S*output;
            normalize_vector(output);

            // sanity check

            assert(form_to_mat(f,output).rank() == 4);

            if(! (output == v)){
                return output;
            }
        }
    }
}

void _find_rank_4(const Form &f, Vector &v, Vector &Out, size_t attempts){
    for (size_t i = 0; i < attempts; i++)
    {
        Matrix M = form_to_mat(f,v);
        M.EF();
        v = M.random_kernel_vector_EF();

        for (auto &u : enumerate_low_rank_neighbours(f, v))
        {
            normalize_vector(u);
            Out = u;
            return;
        }
    }
}


#define FINDER_THREADS 16
#define FINDER_ATTEMPTS 500

Vector find_rank_4(const Form &f){
    std::vector<std::thread> threads;
    std::vector<Vector> vs(FINDER_THREADS,Vector(N));
    std::vector<Vector> outs(FINDER_THREADS);

    std::cout << "Start looking for a rank-4 vector.\n";

    for (size_t i = 0; i < FINDER_THREADS; i++)
    {
        brute_force_rank(f,6,vs[i]);
    }

    while (true)
    {
        for (size_t i = 0; i < FINDER_THREADS; i++)
        {
            threads.emplace_back(_find_rank_4,std::ref(f), std::ref(vs[i]), std::ref(outs[i]), FINDER_ATTEMPTS);
        }
        
        for (size_t i = 0; i < FINDER_THREADS; i++)
        {
            threads[i].join();
        }

        threads.resize(0);

        for (size_t i = 0; i < FINDER_THREADS; i++)
        {
            if(outs[i].size() > 0){
                return outs[i];
            }
        }
    }
}

Vector F(const Form &f, Vector v){
    normalize_vector(v);

    Matrix Mv = form_to_mat(f,v);
    std::vector<Vector> kernel(Mv.kernel());
    
    assert(kernel.size() == 5);

    Matrix MM(kernel.size()*kernel.size(), kernel.size());
    for(size_t kvi=0; kvi < kernel.size(); kvi ++){
        Matrix M = form_to_mat(f,kernel[kvi]);

        // If we need a speedup we can avoid computing half of the rows, because each row corresponds to a symmetric matrix.
        for (size_t i = 0; i < kernel.size(); i++)
        {
            Vector Mkv = M*kernel[i];
            for (size_t j = 0; j < kernel.size(); j++)
            {
                MData(MM, i*kernel.size() + j, kvi) = kernel[j]*Mkv;
            }
        }
    }

    std::vector<Vector> MMKernel = MM.kernel();
    assert(MMKernel.size() == 2);

    Matrix W(6,N);
    W.zero();

    while(1){
        long long r = random() % Q;

        Vector X = MMKernel[0] + r*MMKernel[1];
        
        Vector rad(N,0);
        for (size_t i = 0; i < kernel.size(); i++)
        {
            rad += X[i] * kernel[i]; 
        }
        
        Matrix Mrad = form_to_mat(f,rad);
        Mrad.EF();
        size_t rank = Mrad.rank_EF();
        if(rank == 6){
            auto Mrad_ker = Mrad.kernel_EF();
            for (size_t i = 0; i < N; i++)
            {
                MData(W,3,i) = Mrad_ker[0][i];
                MData(W,4,i) = Mrad_ker[1][i];
                MData(W,5,i) = Mrad_ker[2][i];
            }

            W.EF();
            if(W.rank_EF() == 4){
                break;
            }
        }
    }

    std::vector<Matrix> WMatrices(4);
    for (size_t i = 0; i < 4; i++)
    {
        WMatrices[i] = form_to_mat(f,W.row(i));
    }

    auto sols = enumerate_low_rank(WMatrices);
    if(sols.size() != 2){
        std::cout << "number of vectors in W:" << sols.size() << std::endl;
        return Vector(0);
    }

    for( Vector &sol : sols){
        Vector Fv(N,0);
        for (size_t i = 0; i < 4; i++)
        {
            Fv += sol[i]*W.row(i);
        }
        normalize_vector(Fv);

        if(!(Fv == v)){
            return Fv;
        }
    }

    return Vector(0);
}

std::vector<Vector> compute_F_chain(const Form &f,const Vector &start, size_t Len){
    size_t len = 1;
    std::vector<Vector> FChain(Len);
    FChain[0] = start;
    while( len < Len){
        FChain[len] = F(f,FChain[len-1]);
        //std::cout << len << " | " << FChain[len] << "\n";
        if(FChain[len].size() == 0){
            break;
        }
        len += 1;
    }
    FChain.resize(len);
    return FChain;
}

// returns position on chain of next fresh value
// collumns of frame form the basis
size_t extract_frame(const Form& f, const std::vector<Vector>& chain, size_t pos, Matrix& frame){
    size_t used = 0;
    frame = Matrix(N,N+1);
    Matrix M(N,N);
    frame.zero();
    M.zero();

    // choose basis
    while(used < N){
        if(pos >= chain.size()){
            return 0;
        }
        M.set_row(used, chain[pos]);

        M.EF();
        if(M.rank_EF() == used+1){
            frame.set_col(used, chain[pos]);
            used += 1;
        }
        pos += 1;
    }

    // finish frame
    while (true)
    {
        if(pos+1 >= chain.size()){
            return 0;
        }
        frame.set_col(N,chain[pos]);
        pos += 1;

        Vector kv = frame.kernel()[0];
        bool no_zeros = true;
        for(size_t i=0; i<N; i++){
            if(kv[i]%Q == 0){
                no_zeros = false;
                break;
            } 
            frame.set_col( i, kv[i]*frame.col(i) );
        }
        if(no_zeros){
            break;
        }
    }

    return pos;
}

Vector extract_signature(const Form &f, const std::vector<Vector> &chain, size_t pos){
    Matrix frame;
    pos = extract_frame(f, chain, pos, frame);
    frame.set_col(N,chain[pos]);
    Vector signature = frame.kernel()[0];
    signature.resize(N);
    normalize_vector(signature);
    return signature;
}

Form act_on_form(const Form &f, const Matrix &S){
    //might be off by a factor -1?
    Form Out;
    size_t ctr = 0;
    for (size_t i = 0; i < N; i++)
    {
        Matrix Mi = form_to_mat(f,S.col(i));
        for (size_t j = i+1; j < N; j++)
        {
            Vector MiSj = Mi*S.col(j);
            for (size_t k = j+1; k < N; k++)
            {
                Out[ctr++] =  S.col(k) * MiSj;
            }
        }
    }
    return Out;
}


void worker(std::vector<Vector> &Chain, const Form &f, Vector &start, size_t len){
    start = random_neighbour(f,start);
    start = random_neighbour(f,start);
    start = random_neighbour(f,start);
    Chain = compute_F_chain(f, start, len);
}

#define THREADS 16
#define FCHAIN_LEN1 260
#define FCHAIN_LEN2 30
std::pair<Vector,Vector> find_colliding_pair(const Form &f1, const Form &f2){
    std::unordered_map<Vector,Vector> Sig_map1;
    std::unordered_map<Vector,Vector> Sig_map2;

    long long start_time = time(NULL);
    long long found_rank_4_time;

    Vector v1,v2;
    
    v1 = find_rank_4(f1);
    std::cout << "found first rank-4 vectors after " << time(NULL) - start_time << " seconds.\n";

    v2 = find_rank_4(f2);
    found_rank_4_time = time(NULL);
    std::cout << "found both rank-4 vectors after " << found_rank_4_time - start_time << " seconds.\n";

    std::cout << v1 << std::endl;
    std::cout << v2 << std::endl;

    std::vector<Vector> v1s(THREADS,v1);
    std::vector<Vector> v2s(THREADS,v2);

    while(true)
    {
        {
            // compute signatures in f1
            std::vector<Vector> FChain[THREADS];
            std::vector<std::thread> threads;
            for (size_t i = 0; i < THREADS; i++)
            {
                threads.emplace_back(worker, std::ref(FChain[i]), std::ref(f1), std::ref(v1s[i]), FCHAIN_LEN1 );
            }
            
            for (size_t i = 0; i < THREADS; i++)
            {
                threads[i].join();
            }
            
            for (size_t i = 0; i < THREADS; i++)
            {
                for (size_t pos = 0; pos < FChain[i].size()-N -1; pos++)
                {
                    Vector sig = extract_signature(f1, FChain[i], pos);
                    if(sig.size() == 0){
                        break;
                    }
                    Sig_map1[sig] = FChain[i][pos];
                    if(Sig_map2.count(sig) == 1){
                        std::cout << "found both rank-4 vectors after " << found_rank_4_time - start_time << " seconds.\n";
                        std::cout << "signatures in map1: " << Sig_map1.size() << "\n";
                        std::cout << "signatures in map2: " << Sig_map2.size() << "\n" << std::endl;
                        return std::pair<Vector,Vector>(Sig_map1[sig], Sig_map2[sig]);
                    }
                }
            }
            std::cout << "signatures in map1: " << Sig_map1.size() << "\n";
            std::cout << "signatures in map2: " << Sig_map2.size() << "\n" << std::endl;
        }

        // compute signatures in f2
        for (size_t tr = 0; tr < 8; tr++)
        {
            std::vector<Vector> FChain[THREADS];
            std::vector<std::thread> threads;
            for (size_t i = 0; i < THREADS; i++)
            {
                threads.emplace_back(worker, std::ref(FChain[i]), std::ref(f2), std::ref(v2s[i]), FCHAIN_LEN2 );
            }
        
            for (size_t i = 0; i < THREADS; i++)
            {
                threads[i].join();
            }

            for (size_t i = 0; i < THREADS; i++)
            {            
                Vector sig(0);
                size_t pos;
                for (pos = FChain[i].size()-N-2; pos >= 0; pos--)
                {
                    sig = extract_signature(f2,FChain[i], pos);
                    if(sig.size() > 0){
                        break;
                    }
                }
                
                if(sig.size() > 0){
                    Sig_map2[sig] = FChain[i][pos];
                    if(Sig_map1.count(sig) == 1){
                        std::cout << "found both rank-4 vectors after " << found_rank_4_time - start_time << " seconds.\n";
                        std::cout << "signatures in map1: " << Sig_map1.size() << "\n";
                        std::cout << "signatures in map2: " << Sig_map2.size() << "\n" << std::endl;
                        return std::pair<Vector,Vector>(Sig_map1[sig], Sig_map2[sig]);
                    }
                }
            }
        }
        std::cout << "signatures in map1: " << Sig_map1.size() << "\n";
        std::cout << "signatures in map2: " << Sig_map2.size() << "\n" << std::endl;
    }
}

Matrix find_isomorphsism(const Form &f1, const Vector &v1, const Form &f2, const Vector &v2){

    auto FChain1 = compute_F_chain(f1,v1,N+5);
    auto FChain2 = compute_F_chain(f2,v2,N+5);

    Matrix frame1, frame2;
    extract_frame(f1,FChain1,0,frame1);
    extract_frame(f2,FChain2,0,frame2);

    Matrix basis1(N,N);
    Matrix basis2(N,N);

    for (size_t i = 0; i < N; i++)
    {
        basis1.set_row(i, frame1.col(i));
        basis2.set_row(i, frame2.col(i));
    }

    Matrix isoT(N,N);
    Matrix iso(N,N);

    multiply_by_inverse(basis2, basis1, isoT);
    for (size_t i = 0; i < N; i++)
    {
        iso.set_row(i, isoT.col(i));
    }
    iso.reduce();

    return iso;
}