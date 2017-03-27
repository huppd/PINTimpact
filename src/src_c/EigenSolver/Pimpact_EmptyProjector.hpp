#pragma once
#ifndef PIMPACT_EMPTYPROJECTOR_HPP
#define PIMPACT_EMPTYPROJECTOR_HPP



namespace Pimpact {



template<class OperatorT>
class EmptyProjector {

public:

	EmptyProjector() {}
	EmptyProjector( const Teuchos::RCP<const OperatorT>& op ) {}
	void operator()( typename OperatorT::RangeFieldT& rhs ) const {}

}; // end of class EmptyProjector 

} // end of namespace Pimpact


#endif // end of #ifndef PIMPACT_EMPTYPROJECTOR_HPP
