#ifndef PARTICLE_H
#define PARTICLE_H

#include <cassert>
#include <cstdlib>

#include "SAMRAI/pdat/IndexVariable.h"
#include "SAMRAI/pdat/IndexVariable.C"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Timer.h"
#include "SAMRAI/pdat/CellData.C"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexData.C"
#include "SAMRAI/pdat/IndexDataFactory.h"
#include "SAMRAI/pdat/IndexDataFactory.C"

#include "boost/shared_ptr.hpp"
#include <list>

#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/pdat/CellData.C"
#include "SAMRAI/pdat/IndexData.h"
#include "SAMRAI/pdat/IndexData.C"
#include "SAMRAI/pdat/IndexDataFactory.h"
#include "SAMRAI/pdat/IndexDataFactory.C"


using namespace SAMRAI;
using namespace SAMRAI::tbox;

#define SIGN(x) (((x) < 0.0) ? -1 : ((x) > 0.0))

#define REAL_NUMBERS_IN_A_PARTICLE              \
  (4*PARTICLE_NUMBER_OF_STATES + PARTICLE_REAL_PROPERTIES)
#define INT_NUMBERS_IN_A_PARTICLE 1

#define PARTICLE_NUMBER_OF_STATES 7

#define PARTICLE_INT_PROPERTIES 1
#define PARTICLE_REAL_PROPERTIES 1

class RKParticle
{
 public:
  //  int lv;
  double x_p[PARTICLE_NUMBER_OF_STATES];

  double x_a[PARTICLE_NUMBER_OF_STATES];

  double x_c[PARTICLE_NUMBER_OF_STATES];

  double x_f[PARTICLE_NUMBER_OF_STATES];

  // storing redshift etc.
  double rp[PARTICLE_REAL_PROPERTIES];

  // storing id for particles to
  // tell which is which 
  int ip[PARTICLE_INT_PROPERTIES];
  
  RKParticle()
  {
  }
  ~RKParticle()
  {
  }
  RKParticle(double x_in[], int ip_in[])
  {
    for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      x_p[i] = x_a[i] = x_in[i];
    
    for(int i = 0; i < PARTICLE_INT_PROPERTIES; i++)
      ip[i] = ip_in[i];

    for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i++)
      rp[i] = 0;

  }

  RKParticle(double x_in[], double rp_in[], int ip_in[])
  {
    for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
      x_p[i] = x_a[i] = x_in[i];

    for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i++)
      rp[i] = rp_in[i];
    
    for(int i = 0; i < PARTICLE_INT_PROPERTIES; i++)
      ip[i] = ip_in[i];
      
  }
};

class ParticleContainer
{

public:
 ParticleContainer():
  idx(-9999,-9999,-9999)
   {
   }

   ~ParticleContainer()
   {
     p_list.clear();
   }


   // initializing empty ParticleContainer 
 ParticleContainer(const hier::Index idx_in):
     idx(idx_in)
   {
     
   }

   void addParticle(RKParticle & p)
   {
     p_list.push_back(p);
   }
   
   void copySourceItem(
      const hier::Index& idx_in,
      const hier::IntVector& src_offset,
      const ParticleContainer& src_item)
   {
      NULL_USE(idx);
      NULL_USE(src_offset);
      idx = src_item.idx;
      p_list = src_item.p_list;
   }

   ParticleContainer & operator = (const ParticleContainer & p)
     {
       idx = p.idx;
       p_list = p.p_list;
       return *this;
     }

/* Return the total size that 
 *
 *
 *
 */
   size_t getDataStreamSize()
   {
     // starting with list length
     size_t byte_size = tbox::MessageStream::getSizeof<int>(1);

     // for each particle, 6 double variables needed to store
     // their states.
     byte_size += tbox::MessageStream::getSizeof<double>(REAL_NUMBERS_IN_A_PARTICLE)
       * p_list.size();

     byte_size += tbox::MessageStream::getSizeof<int>(INT_NUMBERS_IN_A_PARTICLE)
       * p_list.size();

     // also adding size for 3d indices
     byte_size += tbox::MessageStream::getSizeof<int>(3);

     
     return byte_size;
   }

/* Pack data into tbox::MessageStream by using stream.pack() function, 
 * see SAMRAI's reference for more detail about tbox::MessageStream.
 *
 *
 */
   void packStream(
      tbox::MessageStream& stream)
   {
     //     const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
     /* if (mpi.getRank() == 0) { */
     /*   TBOX_ERROR("Statistic::packStream error...\n" */
     /*              << "    Processor zero should not pack stat data" << std::endl); */
     /* } */

     stream.pack(&idx[0], 3);
     
     int p_num = p_list.size();
     stream.pack(&p_num, 1);
     
     double temp[PARTICLE_NUMBER_OF_STATES+PARTICLE_INT_PROPERTIES+PARTICLE_REAL_PROPERTIES];

     for (std::list<RKParticle>::iterator it=p_list.begin(); it != p_list.end(); ++it)
     {
       for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i++)
         temp[i] = (*it).x_a[i];

       for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i++)
         temp[PARTICLE_NUMBER_OF_STATES+i] = (*it).rp[i];

       
       for(int i = 0; i < PARTICLE_INT_PROPERTIES; i++)
         temp[PARTICLE_NUMBER_OF_STATES+PARTICLE_REAL_PROPERTIES+i] = (*it).ip[i];

       
       stream.pack(&temp[0], PARTICLE_NUMBER_OF_STATES+PARTICLE_INT_PROPERTIES+PARTICLE_REAL_PROPERTIES);
     }
   }

   void unpackStream(
      tbox::MessageStream& stream,
      const hier::IntVector offset)
   {
     //     const SAMRAI_MPI& mpi(SAMRAI_MPI::getSAMRAIWorld());
     /* if (mpi.getRank() != 0) { */
     /*   TBOX_ERROR("Statistic::unpackStream error...\n" */
     /*              << "    Only processor zero should unpack stat data" << std::endl); */
     /* } */

     stream.unpack(&idx[0], 3);
     
     int incoming_p_num;
     //getting info of number of incoming particles

     stream>>incoming_p_num;

     
     for(int i = 0; i < incoming_p_num; i++)
     {
       double temp_p_info[PARTICLE_NUMBER_OF_STATES+PARTICLE_INT_PROPERTIES+PARTICLE_REAL_PROPERTIES];

       RKParticle temp_p;

       //getting level info
       //       stream>>temp_p_lv_info;

       // getting location and velocity info
       stream.unpack(&temp_p_info[0], PARTICLE_NUMBER_OF_STATES+PARTICLE_INT_PROPERTIES+PARTICLE_REAL_PROPERTIES);

       // copying data to a temp particle object
       for(int i = 0; i < PARTICLE_NUMBER_OF_STATES; i ++)
         temp_p.x_a[i] = temp_p.x_p[i] = temp_p_info[i];
       for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i ++)
         temp_p.rp[i] = temp_p_info[PARTICLE_NUMBER_OF_STATES+i];
       for(int i = 0; i < PARTICLE_INT_PROPERTIES; i ++)
         temp_p.ip[i] = temp_p_info[PARTICLE_NUMBER_OF_STATES+PARTICLE_REAL_PROPERTIES+i];
       // adding this temp particle to the list
       p_list.push_back(temp_p);
     }
   }

   void putToRestart(
      std::shared_ptr<tbox::Database> dbase)
   {
      NULL_USE(dbase);
   }
   void getFromRestart(
      std::shared_ptr<tbox::Database> dbase)
   {
      NULL_USE(dbase);
   }

   // print info in this list
   void print()
   {
     //     int p_num = p_list.size();
     for (std::list<RKParticle>::iterator it=p_list.begin(); it != p_list.end(); ++it)
     {
       
       std::cout<<"("<<(*it).x_a[0]<<", "<<(*it).x_a[1]<<", "<<(*it).x_a[2]<<")  ";
       std::cout<<"("<<(*it).x_a[3]<<", "<<(*it).x_a[4]<<", "<<(*it).x_a[5]<<")\n";
#if EVOLVE_LAMBDA
       std::cout<<(*it).x_a[6]<<" ";
#endif
       std::cout<<sqrt(((*it).x_a[0] - 30) * ((*it).x_a[0] - 30)
                       + ((*it).x_a[1] - 30) * ((*it).x_a[1] - 30)
                       + ((*it).x_a[2] - 30) * ((*it).x_a[2] - 30))<<" ";
       for(int i = 0; i < PARTICLE_REAL_PROPERTIES; i ++)
         std::cout<<(*it).rp[i]<<" ";
       std::cout<<"\n";
     }
     std::cout<<" ";
   }

   /* ParticleContainer& operator=(const ParticleContainer& other) // copy assignment */
   /*   { */
   /*     //       tbox::pout<<idx<<" "<<p_list.size()<<" "<<other.p_list.size()<<"\n"; */
   /*     idx = other.idx; */
   /*     p_list = other.p_list; */
   /*     //       p_list.insert(p_list.end(),other.p_list.begin(), other.p_list.end()); */
   /*   } */
   
   hier::Index idx;   
   std::list<RKParticle> p_list;
};

#endif
