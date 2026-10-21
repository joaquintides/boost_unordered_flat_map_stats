/* Statistical properties of boost::unordered_flat_map and
 * absl::flat_hash_map obtained via simulation.
 * 
 * Expanded to match the theoretical analysis in
 * Cyril Nicaud and Pablo Rotondo, “Fine-Grained Analysis of SIMD-Based Hash
 * Table Implementations.” arXiv preprint arXiv:2610.02385, 2026.
 * https://arxiv.org/abs/2610.02385
 * template<bool RandomProbing=false,bool RandomOverflowBit=false>
 * struct boost_map;
 *   - RandomProbing: probe groups at random (based on hash) instead of
 *     quadratically.
 *   - RandomOverflowBit: the associated overflow bit changes randomly
 *     (based on hash) at each visited group.
 *
 * Copyright 2022-2026 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <iostream>
#include <optional>
#include <random>
#include <vector>
#include <type_traits>

struct pow2_quadratic_prober
{
  pow2_quadratic_prober(std::size_t pos_,std::size_t=0 /* hash */):pos{pos_}{}

  inline std::size_t get()const{return pos;}

  inline bool next(std::size_t mask)
  {
    step+=1;
    pos=(pos+step)&mask;
    return step<=mask;
  }

private:
  std::size_t pos,step=0;
};

struct random_prober
{
  random_prober(std::size_t pos_,std::size_t hash):
    pos{pos_},gen((std::uint_fast32_t)hash){}

  inline std::size_t get()const{return pos;}

  inline bool next(std::size_t mask)
  {
    step+=1;
    pos=std::uniform_int_distribution<std::size_t>(0,mask)(gen);
    return step<=mask;
  }

private:
  std::size_t  pos,step=0;
  std::mt19937 gen;
};

struct find_info
{
  find_info& operator +=(const find_info& x)
  {
    num_hops+=x.num_hops;
    num_cmps+=x.num_cmps;
    return *this;
  }

  std::size_t num_hops=0,num_cmps=0;
};

struct default_overflow_mask_producer
{
  default_overflow_mask_producer(std::size_t hash):
    res(static_cast<unsigned char>(1<<(hash%8))){}

  unsigned char operator()()const{return res;}

private:
  unsigned char res;
};

struct random_overflow_mask_producer
{
  random_overflow_mask_producer(std::size_t hash):
    gen((std::uint_fast32_t)hash){}

  unsigned char operator()()
  {
    return 1<<std::uniform_int_distribution<int>(0,7)(gen);
  }

private:
  std::mt19937 gen;
};

template<bool RandomProbing=false,bool RandomOverflowBit=false>
struct boost_map
{
  static constexpr int N=15;
  using prober=std::conditional_t<
    RandomProbing,random_prober,pow2_quadratic_prober>;
  using overflow_mask_producer=std::conditional_t<
    RandomOverflowBit,
    random_overflow_mask_producer,default_overflow_mask_producer>;

  boost_map(std::size_t capacity):groups{capacity}{}

  bool insert(std::size_t hash)
  {
    overflow_mask_producer omp(hash);
    if(find(hash).second)return false;
    for(prober pb(position_for(hash),hash);;
        pb.next(groups.size()-1)){
      auto& g=groups[pb.get()];
      if(g.n<N){
        g.elements[g.n++]=hash;
        return true;
      }
      else{
        g.overflow|=omp();
      }
    }
  }

  std::pair<find_info,bool> find(std::size_t hash)const
  {
    find_info              info;
    overflow_mask_producer omp(hash);
    for(prober pb(position_for(hash),hash);;
        ++info.num_hops,pb.next(groups.size()-1)){
      auto& g=groups[pb.get()];
      auto first=g.elements.begin(),last=first+g.n;
      auto it=std::find(first,last,hash);
      auto num_cmps=
      info.num_cmps+=std::count_if(first,it,[rh=reduced_hash(hash)](auto x){
        return reduced_hash(x)==rh;
      });
      if(it!=last){
        ++info.num_cmps;
        return {info,true};
      }
      if(!(g.overflow&omp())){
        return {info,false};
      }
    }
  }

  float pr_group_full()const
  {
    auto num_full=std::count_if(
      groups.begin(),groups.end(),[](const auto& g){return g.n==N;});

    return float(num_full)/groups.size();
  }

private:
  static unsigned char reduced_hash(std::size_t hash)
  {
    unsigned char h=(unsigned char)hash;
    return h==0?8:h==1?9:h;
  }

  std::size_t position_for(std::size_t hash)const
  {
    return hash>>shift;
  }

  struct group
  {
    std::array<std::size_t,N> elements;
    std::size_t               n=0;
    unsigned char             overflow=0;
  };

  std::vector<group> groups;
  std::size_t        shift=sizeof(std::size_t)*CHAR_BIT-
                       (std::size_t)(std::bit_width(groups.size()-1));
};

struct abseil_map
{
  static constexpr int N=16;

  abseil_map(std::size_t capacity):elements(capacity*N,std::nullopt){}

  bool insert(std::size_t hash)
  {
    if(find(hash).second)return false;
    std::size_t pos=position_for(hash),
                off=pos%N;
    for(pow2_quadratic_prober pb(pos/N);;pb.next(elements.size()/N-1)){
      auto pos0=pb.get()*N+off;
      auto pos1=find(pos0,pos0+N,std::nullopt);
      if(pos1!=pos0+N){
        element_at(pos1)=hash;
        return true;
      }
    }
  }

  std::pair<find_info,bool> find(std::size_t hash)const
  {
    find_info   info;
    std::size_t pos=position_for(hash),
                off=pos%N;
    for(pow2_quadratic_prober pb(pos/N);;
        ++info.num_hops,pb.next(elements.size()/N-1)){
      auto pos0=pb.get()*N+off,
           pos1=find(pos0,pos0+N,hash);
      info.num_cmps+=num_cmps(pos0,pos1,hash);
      if(pos1!=pos0+N){
        ++info.num_cmps;
        return {info,true};
      }
      if(find(pos0,pos0+N,std::nullopt)!=pos0+N){
        return {info,false};
      }
    }
  }

  float pr_group_full()const
  {
    std::size_t num_full=0;
    for(std::size_t pos=0;pos<elements.size();++pos){
      if(find(pos,pos+N,std::nullopt)==pos+N)++num_full;
    }

    return float(num_full)/elements.size();
  }

private:
  std::size_t position_for(std::size_t hash)const
  {
    return (hash>>7)&(elements.size()-1);
  }

  std::size_t find(
    std::size_t first,std::size_t last,std::optional<std::size_t> x)const
  {
    for(;first!=last;++first){
      if(element_at(first)==x)break;
    }
    return first;
  }

  std::size_t num_cmps(
    std::size_t first,std::size_t last,std::size_t hash)const
  {
    std::size_t res=0;
    for(;first!=last;++first){
      const auto& x=element_at(first);
      if(x && (*x&0x7fu)==(hash&0x7fu))++res;
    }
    return res;
  }

  std::optional<std::size_t>& element_at(std::size_t pos)
  {return elements[pos&(elements.size()-1)];}
  const std::optional<std::size_t>& element_at(std::size_t pos)const
  {return elements[pos&(elements.size()-1)];}

  std::vector<std::optional<std::size_t>> elements;
};

template<typename Map>
void stats_row(std::size_t capacity,float lf)
{
  static constexpr auto N=Map::N;

  Map                                        m(capacity);
  std::size_t                                size=std::size_t(capacity*lf*N);
  std::uniform_int_distribution<std::size_t> dist;
  std::mt19937                               gen(0);
  
  for(std::size_t n=0;n!=size;){
    if(m.insert(dist(gen)))++n;
  }
  
  std::cout<<lf<<";"<<m.pr_group_full()<<";";

  find_info info;
  gen.seed(0);
  for(std::size_t n=0;n!=size;++n){
    info+=m.find(dist(gen)).first;
  }
  std::cout
    <<(!size?0.f:float(info.num_hops)/size)<<";"
    <<(!size?0.f:float(info.num_cmps)/size)<<";";

  info={0,0};
  gen.seed(1);
  for(std::size_t n=0;n!=size;){
    auto [i,b]=m.find(dist(gen));
    if(!b){
      info+=i;
      ++n;
    }
  }
  std::cout
    <<(!size?0.f:float(info.num_hops)/size)<<";"
    <<(!size?0.f:float(info.num_cmps)/size)<<"\n";
}

template<typename Map>
void stats(const char* label)
{
  std::cout
    <<label<<"\n"
    <<"load factor;Pr(group full);"
    <<"E(num hops), successful lookup;E(num cmps), successful lookup;"
    <<"E(num hops), unsuccessful lookup;E(num cmps), unsuccessful lookup\n";

  constexpr int   num_points=101;
  constexpr float mlf=0.875f;

  for(int i=0;i<num_points;++i){
    stats_row<Map>(0x20000ul,mlf*i/(num_points-1));
  }
}

int main()
{
  stats<boost_map<>>("boost::unordered_flat_map");
  stats<boost_map<true>>("boost::unordered_flat_map, random probing");
  stats<boost_map<false,true>>(
    "boost::unordered_flat_map, random overflow bit");
  stats<boost_map<true,true>>(
    "boost::unordered_flat_map, random probing and overflow bit");
  stats<abseil_map>("absl::flat_hash_map");
}
