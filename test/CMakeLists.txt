
#set(CMAKE_INCLUDE_CURRENT_DIR ../src)
include_directories("../src/src_c")

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

ADD_EXECUTABLE(ScalarModeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_ScalarModeField_test.cpp )
TARGET_LINK_LIBRARIES( ScalarModeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME ScalarModeFieldTest COMMAND mpirun -np 4 ScalarModeFieldTest)

ADD_EXECUTABLE(VectorModeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_VectorModeField_test.cpp )
TARGET_LINK_LIBRARIES( VectorModeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME VectorModeFieldTest COMMAND mpirun -np 4 VectorModeFieldTest)

ADD_EXECUTABLE(CompoundFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_CompoundField_test.cpp )
TARGET_LINK_LIBRARIES( CompoundFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME CompoundFieldTest COMMAND mpirun -np 4 CompoundFieldTest)

ADD_EXECUTABLE(MultiFieldScalarTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_Scalar_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldScalarTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME MultiFieldScalarTest COMMAND mpirun -np 4 MultiFieldScalarTest)

ADD_EXECUTABLE(MultiFieldVectorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_Vector_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldVectorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME MultiFieldVectorTest COMMAND mpirun -np 4 MultiFieldVectorTest)

ADD_EXECUTABLE(MultiFieldScalarModeTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_ScalarMode_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldScalarModeTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME MultiFieldScalarModeTest COMMAND mpirun -np 4 MultiFieldScalarModeTest)

ADD_EXECUTABLE(MultiFieldVectorModeTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_VectorMode_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldVectorModeTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME MultiFieldVectorModeTest COMMAND mpirun -np 4 MultiFieldVectorModeTest)

ADD_EXECUTABLE(BelosMVOPTest Teuchos_StandardParallelUnitTestMain.cpp BelosPimpact_MVOP_test.cpp )
TARGET_LINK_LIBRARIES( BelosMVOPTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosMVOPTest COMMAND mpirun -np 4 BelosMVOPTest)

ADD_EXECUTABLE(OperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Operator_test.cpp )
TARGET_LINK_LIBRARIES( OperatorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME OperatorTest COMMAND mpirun -np 4 OperatorTest)

ADD_EXECUTABLE(BelosSolverTest Teuchos_StandardParallelUnitTestMain.cpp BelosSolvesPimpact_test.cpp )
TARGET_LINK_LIBRARIES( BelosSolverTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosSolverTest COMMAND mpirun -np 4 BelosSolverTest)

ADD_EXECUTABLE(BelosSolverTestBC Teuchos_StandardParallelUnitTestMain.cpp BelosSolvesPimpactBC_test.cpp )
TARGET_LINK_LIBRARIES( BelosSolverTestBC ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosSolverTestBC COMMAND mpirun -np 4 BelosSolverTestBC)

ADD_EXECUTABLE(LinearProblem Teuchos_StandardParallelUnitTestMain.cpp Pimpact_LinearProblem_test.cpp )
TARGET_LINK_LIBRARIES( LinearProblem ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME LinearProblem COMMAND mpirun -np 4 LinearProblem)

ADD_EXECUTABLE( NOXVectorTest Teuchos_StandardParallelUnitTestMain.cpp NOX_Pimpact_Vector_test.cpp )
TARGET_LINK_LIBRARIES( NOXVectorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NOXVectorTest COMMAND mpirun -np 4 NOXVectorTest)
#ADD_Test(mytets MyUnitTests)

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
