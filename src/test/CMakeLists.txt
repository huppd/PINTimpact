include_directories("../src/src_c")
include_directories("../src/src_c/Space")
include_directories("../src/src_c/Grid")
include_directories("../src/src_c/Fields")
include_directories("../src/src_c/Operators")
include_directories("../src/src_c/Operators/BasicOperators")
include_directories("../src/src_c/Operators/CompoundOperators")
include_directories("../src/src_c/Operators/ModeOperators")
include_directories("../src/src_c/Operators/MultiHarmonicOperators")
include_directories("../src/src_c/Operators/TimeOperators")
include_directories("../src/src_c/BelosPimpact")
include_directories("../src/src_c/NOXPimpact")
include_directories("../src/src_c/MultiGrid")

enable_testing()

ADD_EXECUTABLE( SpaceGridTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_SpaceGrid_test.cpp )
TARGET_LINK_LIBRARIES( SpaceGridTest  ${Pimpact_TPL_LIBRARIES} )
add_test(NAME SpaceGridTest COMMAND mpirun -np 4 SpaceGridTest)

ADD_EXECUTABLE(ScalarFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_ScalarField_test.cpp )
TARGET_LINK_LIBRARIES(ScalarFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE(VectorFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_VectorField_test.cpp )
TARGET_LINK_LIBRARIES( VectorFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( MultiFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( MultiHarmonicFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiHarmonicField_test.cpp )
TARGET_LINK_LIBRARIES( MultiHarmonicFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE(OperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Operator_test.cpp )
TARGET_LINK_LIBRARIES( OperatorTest ${Pimpact_TPL_LIBRARIES} ) 

foreach(dim RANGE 2 3)
	foreach(domain RANGE 1 5)
		add_test(NAME "ScalarFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 8 ScalarFieldTest --domain=${domain} --dim=${dim} )
		add_test(NAME "VectorFieldTest_domain_domain${domain}_dim${dim}" COMMAND mpirun -np 8 VectorFieldTest --domain=${domain} --dim=${dim} ) 
		add_test( NAME "MultiFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 4 MultiFieldTest --domain=${domain} --dim=${dim} )
		add_test(NAME "MultiHarmonicFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 4 MultiHarmonicFieldTest --domain=${domain} --dim=${dim} )
	endforeach( domain )
	add_test(NAME "OperatorTest_dim${dim}" COMMAND mpirun -np 4 OperatorTest --dim=${dim} )
endforeach( dim )

ADD_EXECUTABLE( TimeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_TimeField_test.cpp )
TARGET_LINK_LIBRARIES( TimeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME TimeFieldTest COMMAND mpirun -np 8 TimeFieldTest )

#ADD_EXECUTABLE( ModeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_ModeField_test.cpp )
#TARGET_LINK_LIBRARIES( ModeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
#add_test(NAME ModeFieldTest COMMAND mpirun -np 4 ModeFieldTest)

#ADD_EXECUTABLE(CompoundFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_CompoundField_test.cpp )
#TARGET_LINK_LIBRARIES( CompoundFieldTest ${Pimpact_TPL_LIBRARIES} ) 
#add_test(NAME CompoundFieldTest COMMAND mpirun -np 4 CompoundFieldTest)

ADD_EXECUTABLE(BelosMVOPTest Teuchos_StandardParallelUnitTestMain.cpp BelosPimpact_MVOP_test.cpp )
TARGET_LINK_LIBRARIES( BelosMVOPTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME BelosMVOPTest COMMAND mpirun -np 4 BelosMVOPTest)


ADD_EXECUTABLE(NonlinearOperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_NonlinearOperator_test.cpp )
TARGET_LINK_LIBRARIES( NonlinearOperatorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test(NAME NonlinearOperatorTest COMMAND mpirun -np 4 NonlinearOperatorTest)

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

ADD_EXECUTABLE( MultiGridTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiGrid_test.cpp )
TARGET_LINK_LIBRARIES( MultiGridTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME MultiGridTest COMMAND mpirun -np 4 MultiGridTest )