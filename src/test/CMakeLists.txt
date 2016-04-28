include_directories( "../src_c" )
include_directories( "../src_c/Space" )
include_directories( "../src_c/Fields" )
include_directories( "../src_c/Operators" )
include_directories( "../src_c/Operators/BasicOperators" )
include_directories( "../src_c/Operators/CompoundOperators" )
include_directories( "../src_c/Operators/ModeOperators" )
include_directories( "../src_c/Operators/MultiHarmonicOperators" )
include_directories( "../src_c/Operators/TimeOperators" )
include_directories( "../src_c/EigenSolver" )
include_directories( "../src_c/BelosPimpact" )
include_directories( "../src_c/NOXPimpact" )
include_directories( "../src_c/MultiGrid" )
include_directories( "../src_c/MultiGrid/Transfers" )

enable_testing()

ADD_EXECUTABLE( SpaceGridTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_SpaceGrid_test.cpp )
TARGET_LINK_LIBRARIES( SpaceGridTest  ${Pimpact_TPL_LIBRARIES} )
add_test( NAME SpaceGridTest COMMAND mpirun -np 4 ./SpaceGridTest )

ADD_EXECUTABLE( FieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Field_test.cpp )
TARGET_LINK_LIBRARIES( FieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( VectorFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_VectorField_test.cpp )
TARGET_LINK_LIBRARIES( VectorFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( MultiFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiField_test.cpp )
TARGET_LINK_LIBRARIES( MultiFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( MultiHarmonicFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiHarmonicField_test.cpp )
TARGET_LINK_LIBRARIES( MultiHarmonicFieldTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( OperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_Operator_test.cpp )
TARGET_LINK_LIBRARIES( OperatorTest ${Pimpact_TPL_LIBRARIES} ) 

ADD_EXECUTABLE( EigenSolverTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_EigenSolver_test.cpp )
TARGET_LINK_LIBRARIES( EigenSolverTest ${Pimpact_TPL_LIBRARIES} ) 

foreach(dim RANGE 3 3)
	foreach( domain RANGE 0 4 )
		add_test( NAME "FieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 8 ./FieldTest --domain=${domain} --dim=${dim} )
		add_test( NAME "VectorFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 8 ./VectorFieldTest --domain=${domain} --dim=${dim} ) 
		add_test( NAME "MultiFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 8 ./MultiFieldTest --domain=${domain} --dim=${dim} )
		add_test( NAME "MultiHarmonicFieldTest_domain${domain}_dim${dim}" COMMAND mpirun -np 8 ./MultiHarmonicFieldTest --domain=${domain} --dim=${dim} )
	endforeach( domain )
	add_test( NAME "OperatorTest_dim${dim}" COMMAND mpirun -np 8 ./OperatorTest --dim=${dim} --npx=2 --npy=2 --npz=2 )
	add_test( NAME "EigenSolverTest_dim${dim}" COMMAND mpirun -np 8 ./EigenSolverTest --dim=${dim} --npx=2 --npy=2 --npz=2 )
endforeach( dim )


ADD_EXECUTABLE( TimeFieldTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_TimeField_test.cpp )
TARGET_LINK_LIBRARIES( TimeFieldTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME TimeFieldTest COMMAND mpirun -np 8 ./TimeFieldTest )

ADD_EXECUTABLE( SmootherTest Teuchos_StandardParallelUnitTestMain.cpp
Pimpact_Smoother_test.cpp )
TARGET_LINK_LIBRARIES( SmootherTest ${Pimpact_TPL_LIBRARIES} )
add_test( NAME SmootherTest COMMAND mpirun -np 8 SmootherTest )

ADD_EXECUTABLE( BelosMVOPTest Teuchos_StandardParallelUnitTestMain.cpp BelosPimpact_MVOP_test.cpp )
TARGET_LINK_LIBRARIES( BelosMVOPTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME BelosMVOPTest COMMAND mpirun -np 4 ./BelosMVOPTest )


ADD_EXECUTABLE( NonlinearOperatorTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_NonlinearOperator_test.cpp )
TARGET_LINK_LIBRARIES( NonlinearOperatorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME NonlinearOperatorTest COMMAND mpirun -np 8 ./NonlinearOperatorTest )

ADD_EXECUTABLE( LinearProblemTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_LinearProblem_test.cpp )
TARGET_LINK_LIBRARIES( LinearProblemTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME LinearProblemTest COMMAND mpirun -np 4 ./LinearProblemTest )

ADD_EXECUTABLE( NOXVectorTest Teuchos_StandardParallelUnitTestMain.cpp NOX_Pimpact_Vector_test.cpp )
TARGET_LINK_LIBRARIES( NOXVectorTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME NOXVectorTest COMMAND mpirun -np 4 ./NOXVectorTest )

ADD_EXECUTABLE( NOXSolvesPimpactTest Teuchos_StandardParallelUnitTestMain.cpp NOXSolvesPimpact_test.cpp )
TARGET_LINK_LIBRARIES( NOXSolvesPimpactTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME NOXSolvesPimpactTest COMMAND mpirun -np 4 ./NOXSolvesPimpactTest )

ADD_EXECUTABLE( MultiGridTest Teuchos_StandardParallelUnitTestMain.cpp Pimpact_MultiGrid_test.cpp )
TARGET_LINK_LIBRARIES( MultiGridTest ${Pimpact_TPL_LIBRARIES} ) 
add_test( NAME MultiGridTest COMMAND mpirun -np 8 MultiGridTest )

ADD_EXECUTABLE( MultiGridTimeTest Teuchos_StandardParallelUnitTestMain.cpp
Pimpact_MultiGridTime_test.cpp )
TARGET_LINK_LIBRARIES( MultiGridTimeTest ${Pimpact_TPL_LIBRARIES} )
add_test( NAME MultiGridTimeTest COMMAND mpirun -np 8 MultiGridTimeTest )
