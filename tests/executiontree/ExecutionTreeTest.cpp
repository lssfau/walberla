#include "executiontree/ExecutionTree.h"

#include <iostream>
#include "core/logging/Logging.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

using namespace walberla;
namespace et = executiontree;

class MyFunctor
{
public:
   void operator() ()
   {
      WALBERLA_LOG_RESULT( "i = " << i );
      i += 1;
   }

   int i = 0;
};


int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   int counter1 = 0;
   auto func1 = [&counter1]() {
      WALBERLA_LOG_RESULT("A");
      ++counter1;
   };

   int counter2 = 0;
   auto func2 = [&counter2]() {
      ++counter2;
   };

   int counter3 = 0;
   auto func3 = [&counter3]() {
      ++counter3;
   };

   auto func4 = [] {  WALBERLA_LOG_RESULT("B"); };

   auto myFunctor = make_shared<MyFunctor>();

   auto s = et::parallelSequence( { et::everyNth( et::functor( func2, "func2" ), 5, true ),
                                    et::everyNth( et::functor( func3, "func3" ), 5, false ),
                                    et::functor( func1, "func1" ),
                                    et::functor( func4, "func4" ),
                                    et::functor( myFunctor, "myFunctor") } );

   auto l = et::loop( s, 20 );
   myFunctor->i = 42;

   std::cout << *l << std::endl;
   l->run();

   WALBERLA_CHECK_EQUAL( counter1, 20 );
   WALBERLA_CHECK_EQUAL( counter2, 20 / 5 );
   WALBERLA_CHECK_EQUAL( counter3, ( 20 / 5 ) - 1 );
   WALBERLA_CHECK_EQUAL( myFunctor->i, 20 + 42 );

   return 0;
}