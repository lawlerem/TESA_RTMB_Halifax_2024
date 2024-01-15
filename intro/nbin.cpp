#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  PARAMETER(logsize);

  PARAMETER(logitp);
  Type size = exp(logsize);
  Type p=invlogit(logitp);
  Type nll = -sum(dnbinom(Y, size, p, true));
  ADREPORT(size);
  return nll;
}

