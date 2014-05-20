
#set(CMAKE_INCLUDE_CURRENT_DIR ../src)
include_directories("../src/src_c")
include_directories("../src/src_c/Space")
include_directories("../src/src_c/Grid")
include_directories("../src/src_c/Fields")
include_directories("../src/src_c/Operators")
include_directories("../src/src_c/NOX_Pimpact")

enable_testing()

ADD_EXECUTABLE( SpaceTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Space_test.cpp )
TARGET_LINK_LIBRARIES( SpaceTest  ${Pimpact_TPL_LIBRARIES} )
add_test(NAME SpaceTest COMMAND mpirun -np 4 SpaceTest)

ADD_EXECUTABLE(ScalarFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_ScalarField_test.cpp )
TARGET_LINK_LIBRARIES(ScalarFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME ScalarFieldTest COMMAND mpirun -np 4 ScalarFieldTest)

ADD_EXECUTABLE(VectorFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_VectorField_test.cpp )
TARGET_LINK_LIBRARIES( VectorFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME VectorFieldTest COMMAND mpirun -np 4 VectorFieldTest)

ADD_EXECUTABLE( ModeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_ModeField_test.cpp )
TARGET_LINK_LIBRARIES( ModeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME ModeFieldTest COMMAND mpirun -np 4 ModeFieldTest)

ADD_EXECUTABLE(CompoundFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_CompoundField_test.cpp )
TARGET_LINK_LIBRARIES( CompoundFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME CompoundFieldTest COMMAND mpirun -np 4 CompoundFieldTest)

ADD_EXECUTABLE( MultiHarmonicFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiHarmonicField_test.cpp )
TARGET_LINK_LIBRARIES( MultiHarmonicFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME MultiHarmonicFieldTest COMMAND mpirun -np 4 MultiHarmonicFieldTest)

ADD_EXECUTABLE( MultiFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME MultiFieldTest COMMAND mpirun -np 4 MultiFieldTest )

ADD_EXECUTABLE(BelosMVOPTest Teuchos_StandardParallelUnitTestMain.cpp BelosPimpact_MVOP_test.cpp )
TARGET_LINK_LIBRARIES( BelosMVOPTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosMVOPTest COMMAND mpirun -np 4 BelosMVOPTest)

ADD_EXECUTABLE(OperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Operator_test.cpp )
TARGET_LINK_LIBRARIES( OperatorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME OperatorTest COMMAND mpirun -np 4 OperatorTest)

ADD_EXECUTABLE(BelosSolverTest Teuchos_StandardParallelUnitTestMain.cpp BelosSolvesPimpact_test.cpp )
TARGET_LINK_LIBRARIES( BelosSolverTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosSolverTest COMMAND mpirun -np 4 BelosSolverTest)

# ADD_EXECUTABLE(BelosSolverTestBC Teuchos_StandardParallelUnitTestMain.cpp BelosSolvesPimpactBC_test.cpp )
# TARGET_LINK_LIBRARIES( BelosSolverTestBC ${Pimpact_TPL_LIBRARIES} ) 
# add_test(NAME BelosSolverTestBC COMMAND mpirun -np 4 BelosSolverTestBC)

ADD_EXECUTABLE(LinearProblemTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_LinearProblem_test.cpp )
TARGET_LINK_LIBRARIES( LinearProblemTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME LinearProblemTest COMMAND mpirun -np 4 LinearProblemTest)

ADD_EXECUTABLE( NOXVectorTest Teuchos_StandardParallelUnitTestMain.cpp NOX_Pimpact_Vector_test.cpp )
TARGET_LINK_LIBRARIES( NOXVectorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NOXVectorTest COMMAND mpirun -np 4 NOXVectorTest)

ADD_EXECUTABLE( NOXInterfaceTest Teuchos_StandardParallelUnitTestMain.cpp NOX_Pimpact_Interface_test.cpp )
TARGET_LINK_LIBRARIES( NOXInterfaceTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NOXInterfaceTest COMMAND mpirun -np 4 NOXInterfaceTest)

ADD_EXECUTABLE( NOXGroupTest Teuchos_StandardParallelUnitTestMain.cpp NOX_Pimpact_Group_test.cpp )
TARGET_LINK_LIBRARIES( NOXGroupTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NOXGroupTest COMMAND mpirun -np 4 NOXGroupTest)

ADD_EXECUTABLE( NOXSolvesPimpactTest Teuchos_StandardParallelUnitTestMain.cpp NOXSolvesPimpact_test.cpp )
TARGET_LINK_LIBRARIES( NOXSolvesPimpactTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NOXSolvesPimpactTest COMMAND mpirun -np 4 NOXSolvesPimpactTest)

#INCLUDE(TribitsAddExecutableAndTest)

#project("")
#Tribits_ADD_EXECUTABLE_AND_TEST(
  #Simple_UnitTests
  #SOURCES
	#Int_UnitTests.cpp
	##TemplateFunc_UnitTests.cpp
	##vector_UnitTests.cpp
	#${TEUCHOS_STD_UNIT_TEST_MAIN}
  #NUM_MPI_PROCS 1
  #STANDARD_PASS_OUTPUT
  #)
#pimpact_ADD_EXECUTABLE_AND_TEST(
#ADD_EXECUTABLE(
#set(pimpact_TEST_CATEGORIES hello)
#Tribits_ADD_EXECUTABLE_AND_TEST(
  #MyUnitTests
  #SOURCES
	#UnitTestMain.cpp
	#Int_UnitTests.cpp
	##... other unit test files ...
  #STANDARD_PASS_OUTPUT
  #)
#ADD_EXECUTABLE(MyUnitTests ${TEUCHOS_STD_UNIT_TEST_MAIN} UnitTestMain.cpp )
#ADD_EXECUTABLE(MyUnitTests ${TEUCHOS_STD_UNIT_TEST_MAIN} UnitTestMain.cpp )
