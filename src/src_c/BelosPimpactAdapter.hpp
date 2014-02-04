//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#pragma once
#ifndef BELOSPIMPACTADAPTER_HPP
#define BELOSPIMPACTADAPTER_HPP


/*! \file BelosPimpactAdapter.hpp
    \brief Provides several interfaces between Belos virtual classes and Pimpact concrete classes.
*/


#include <Pimpact_MultiField.hpp>
#include <Pimpact_OperatorMV.hpp>

//#include <Teuchos_Assert.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
//#include <Teuchos_Array.hpp>
//#include <Teuchos_DefaultSerialComm.hpp>

//#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

#ifdef HAVE_BELOS_TSQR
#	include "BelosStubTsqrAdapter.hpp"
#endif // HAVE_BELOS_TSQR


namespace Belos {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of \c Belos::MultiVecTraits for \c Pimpact::MultiField.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of MultiVecTraits for \c MV = \c Pimpact::MultiField.
  ///
  /// This interface lets Belos' solvers work directly with
  ///
  /// The two template parameters of this partial specialization
  /// correspond to the \c Scalar type, and to the inner \c Field type of \c Pimpact::MultiField
  /// \c Pimpact::MultiField.  See the \c Pimpact::MultiField documentation
  /// for more information.
  template<class Scalar, class Field>
  class MultiVecTraits<Scalar, Pimpact::MultiField<Field> > {
    typedef Pimpact::MultiField<Field> MV;
  public:
//#ifdef HAVE_BELOS_Pimpact_TIMERS
//    static Teuchos::RCP<Teuchos::Time> mvTimesMatAddMvTimer_, mvTransMvTimer_;
//#endif

    /// \brief Create a new multivector with \c numvecs columns.
    static Teuchos::RCP<Pimpact::MultiField<Field> >
    Clone (const Pimpact::MultiField<Field>& mv, const int numvecs)
    {
    	return mv.Clone(numvecs);
    }


    /// \brief Create a new deep copy of multivector.
    static Teuchos::RCP<Pimpact::MultiField<Field> > CloneCopy( const Pimpact::MultiField<Field>& mv )
    {
    	return mv.CloneCopy();
    }


    /// \brief Create a new deep copy of multivector, considering only the \c Field's of \c index.
    static Teuchos::RCP<Pimpact::MultiField<Field> > CloneCopy( const Pimpact::MultiField<Field>& mv, const std::vector<int>& index )
    {
    	return mv.CloneCopy(index);
    }


    /// \brief Create a new deep copy of multivector, considering only the \c Field's in the range of \c index.
    static Teuchos::RCP<Pimpact::MultiField<Field> >
    CloneCopy (const Pimpact::MultiField<Field>& mv,
               const Teuchos::Range1D& index)
    {
      return mv.CloneCopy(index);
    }


    /// \brief Create a new view to multivector, considering only the \c Field's of \c index.
    static Teuchos::RCP<Pimpact::MultiField<Field> >
    CloneViewNonConst (Pimpact::MultiField<Field>& mv,
                       const std::vector<int>& index)
    {
    	return mv.CloneViewNonConst(index);
    }


    /// \brief Create a new view to multivector, considering only the \c Field's in the range of \c index.
    static Teuchos::RCP<Pimpact::MultiField<Field> >
    CloneViewNonConst (Pimpact::MultiField<Field>& mv,
                       const Teuchos::Range1D& index)
    {
    	return mv.CloneViewNonConst(index);
    }


    /// \brief Create a new const view to multivector, considering only the \c Field's of \c index.
    static Teuchos::RCP<const Pimpact::MultiField<Field> >
    CloneView (const Pimpact::MultiField<Field>& mv,
               const std::vector<int>& index)
    {
    	return mv.CloneView(index);
    }


    /// \brief Create a new const view to multivector, considering only the \c Field's in the range of \c index.
    static Teuchos::RCP<const Pimpact::MultiField<Field> >
    CloneView (const Pimpact::MultiField<Field>& mv,
               const Teuchos::Range1D& index)
    {
    	return mv.CloneView(index);
    }


    /// \brief return the number of the Vector/Field length, it is assumed that every \c Field of the multivector has the same.
    static int GetVecLength( const Pimpact::MultiField<Field>& mv )
    { return mv.GetVecLength(); }


    /// \brief return the number of the Vector/Field's.
    static int GetNumberVecs( const Pimpact::MultiField<Field>& mv )
    { return mv.GetNumberVecs(); }

    static bool HasConstantStride( const Pimpact::MultiField<Field>& mv )
    { return true; }

    /// \brief <tt> mv:= alpha A*B + beta mv </tt>
    static void
    MvTimesMatAddMv (const Scalar& alpha,
                     const Pimpact::MultiField<Field>& A,
                     const Teuchos::SerialDenseMatrix<int,Scalar>& B,
                     const Scalar& beta,
                     Pimpact::MultiField<Field>& mv)
    {
    	mv.TimesMatAdd(alpha,A,B,beta);
    }

    /// \brief <tt>mv := alpha*A + beta*B</tt>
    static void
    MvAddMv (Scalar alpha,
             const Pimpact::MultiField<Field>& A,
             Scalar beta,
             const Pimpact::MultiField<Field>& B,
             Pimpact::MultiField<Field>& mv)
    {
      mv.Add( alpha, A, beta, B );
    }


    /// \brief <tt>mv := alpha*mv </tt>
    static void
    MvScale (Pimpact::MultiField<Field>& mv,
             const Scalar& alpha)
    {
      mv.Scale (alpha);
    }


    /// \brief <tt>mv[i] := alpha[i]*mv[i] </tt>
    static void
    MvScale (Pimpact::MultiField<Field>& mv,
             const std::vector<Scalar>& alphas)
    {
      mv.Scale (alphas);
    }


    /// \brief <tt>C[j,i] := dot(A[j],B[i]) </tt>
    static void
    MvTransMv (Scalar alpha,
               const Pimpact::MultiField<Field>& A,
               const Pimpact::MultiField<Field>& B,
               Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
    	A.Trans(alpha,B,C);
    }


    /// \brief For all columns j of A, set <tt>dots[j] := A[j]^T * B[j]</tt>.
    static void
    MvDot (const Pimpact::MultiField<Field>& A,
           const Pimpact::MultiField<Field>& B,
           std::vector<Scalar>& dots)
    {
    	A.Dot(B,dots);
    }


    /// \brief For all columns j of mv, set <tt>normvec[j] = norm(mv[j])</tt>.
    static void
    MvNorm (const Pimpact::MultiField<Field>& mv,
            std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec,
            NormType type=TwoNorm)
    {
    	mv.Norm(normvec,type);
    }


    /// \brief <tt> mv[i]:=A[index[i]] </tt>
    static void SetBlock( const Pimpact::MultiField<Field>& A,
    		const std::vector<int>& index,
    		Pimpact::MultiField<Field>& mv )
    {
    	mv.SetBlock(A,index);
    }


    /// \brief <tt> mv[i]:=A[i-index.lb] </tt>
    static void
    SetBlock (const Pimpact::MultiField<Field>& A,
              const Teuchos::Range1D& index,
              Pimpact::MultiField<Field>& mv)
    {
    	mv.SetBlock(A,index);
    }


    /// \brief <tt> mv:=A </tt>
    static void
    Assign (const Pimpact::MultiField<Field>& A,
            Pimpact::MultiField<Field>& mv)
    {
    	mv.Assign(A);
    }


    /// \brief make \c mv uniform random (0,1)
    static void MvRandom( Pimpact::MultiField<Field>& mv )
    {
      mv.Random();
    }

    /// \brief init \c mv:= \c alpha everywhere(inner field)
    static void MvInit( Pimpact::MultiField<Field>& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.Init(alpha); }

    /// \brief print function(for debbuging)
    /// \todo add \c os to mv.Print()
    static void MvPrint( const Pimpact::MultiField<Field>& mv, std::ostream& os )
    {
//      Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
      mv.Print();
    }

#ifdef HAVE_BELOS_TSQR
//    /// \typedef tsqr_adaptor_type
//    /// \brief stub TsqrAdaptor specialization for Pimpact::MultiField
//    ///
    typedef Belos::details::StubTsqrAdapter< Pimpact::MultiField<Field> > tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  };

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Pimpact::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of \c Belos::OperatorTraits for \c Pimpact::OperatorMV.
  /// it has three template parameters Scalar, Field, and the inner Operator
  /// \note sadly Belos allows only Operators having the same domain and range
  template <class Scalar, class Field, class Operator>
  class OperatorTraits <Scalar, Pimpact::MultiField<Field>, Pimpact::OperatorMV<Operator> >
  {
  public:
  	/// \brief applys the inner \c Operator, such that \c Y:= \c Op( \c X)
  	/// \note up to now only no NOTRANS operators can be handled
    static void
    Apply (const Pimpact::OperatorMV<Operator>& Op,
           const Pimpact::MultiField<Field>& X,
           Pimpact::MultiField<Field>& Y,
           Belos::ETrans trans=NOTRANS) {
//				std::cout << "x.getVecLength(): " << X.GetVecLength()<< "\n";
//      		std::cout << "y.getVecLength(): " << Y.GetVecLength()<< "\n";
          Op.apply(X,Y,NOTRANS);
    }
//      switch (trans) {
//        case NOTRANS:
//          Op.apply(X,Y,NOTRANS);
//          break;
//        case TRANS:
////          Op.apply(X,Y,TRANS);
//          break;
//        case CONJTRANS:
////          Op.apply(X,Y,CONJTRANS);
//          break;
//      default:
//        const std::string scalarName = Teuchos::TypeNameTraits<Scalar>::name();
//        const std::string loName = Teuchos::TypeNameTraits<LO>::name();
//        const std::string goName = Teuchos::TypeNameTraits<GO>::name();
//        const std::string nodeName = Teuchos::TypeNameTraits<Node>::name();
//        const std::string otName = "Belos::OperatorTraits<" + scalarName
//          + "," + loName + "," + goName + "," + nodeName + ">";
//        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, otName << ": Should never "
//                           "get here; fell through a switch statement.  "
//                           "Please report this bug to the Belos developers.");
//      }

    /**
     * @param Op
     * @return \c true if Op has a transpose implemented
     */
    static bool
    HasApplyTranspose (const Pimpact::OperatorMV<Operator>& Op)
    {
      return Op.hasTransposeApply();
    }
  };
//
//  // Partial specialization for MV=Pimpact::MultiField.
//  template<class Scalar, class LO, class GO, class Node>
//  class MultiVecTraitsExt<Scalar, Pimpact::MultiField<Scalar, LO, GO, Node> > {
//  public:
//    typedef Pimpact::MultiField<Scalar, LO, GO, Node> MV;
//    static ptrdiff_t GetGlobalLength( const MV& mv ) {
//      return Teuchos::as<ptrdiff_t> (mv.getGlobalLength ());
//    }
//  };

} // end of Belos namespace

#endif // end of file BELOSPIMPACTADAPTER_HPP
