
#include "matrix.h"
#include "minrank.h"
#include "form.h"

#include <unordered_map>

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

void find_loops(const Form &f){

    Vector v = find_rank_4(f);

    while(true){

        v = random_neighbour(f,v);
        v = random_neighbour(f,v);
        v = random_neighbour(f,v);

        Vector u = v;
        long long count = 0;
        std::unordered_map<Vector,long long> Map;

        while (true)
        {
            Map[v] = count;
            count ++;
            Vector Fv = F(f,v);
            if(Map.count(Fv) == 1){
                std::cout << "level of starting node: " << Map[Fv] << std::endl;
                break;
            }
        }
    }
}

int main(int argc, char const *argv[])
{

    srand(time(NULL));

    // generate f1, f2
    std::cout << "Generating challenge ... \n";
    Form f1; 
    randomize_form(f1);


    Matrix S = random_invertible_matrix(N);
    Form f2 = act_on_form(f1,S);
    S.zero(); // forget S

    std::cout << "Form 1: ";
    print_array(f1);
    std::cout << "Form 2: ";
    print_array(f2);

    long long start_time = time(NULL);
    // find isomorphism
    auto collision = find_colliding_pair(f1, f2);
    std::cout << "colliding pair found:" << std::endl ;
    std::cout << collision.first << std::endl ;
    std::cout << collision.second << std::endl ;

    std::cout << "It took " << time(NULL)-start_time << " seconds to find a collision.\n";

    Matrix iso = find_isomorphsism(f1, collision.first, f2, collision.second);

    std::cout << "the total attack took: " << time(NULL)-start_time << " seconds.\n\n";

    std::cout << "Isomorphism: \n" << iso ;

    // sanity check
    Form f2_prime = act_on_form(f1, iso);
    std::cout << "f2:       " ;
    print_array(f2);
    std::cout << "f2_prime: " ;
    print_array(f2_prime);

    normalize_array(f2);
    normalize_array(f2_prime);

    std::cout << "normalized:\n";
    std::cout << "f2:       " ;
    print_array(f2);
    std::cout << "f2_prime: " ;
    print_array(f2_prime);

    return 0;
}