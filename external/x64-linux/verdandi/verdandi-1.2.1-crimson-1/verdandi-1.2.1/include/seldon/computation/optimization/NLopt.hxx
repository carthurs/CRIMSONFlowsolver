// Copyright (c) 2007-2010 Massachusetts Institute of Technology
// Modifications by Marc Fragu, copyright (C) 2011 INRIA
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.


#ifndef SELDON_COMPUTATION_OPTIMIZATION_NLOPT_HXX
#define SELDON_COMPUTATION_OPTIMIZATION_NLOPT_HXX


#include <nlopt.hpp>


namespace nlopt
{


  typedef double (*svfunc)(const Seldon::Vector<double> &x,
                           Seldon::Vector<double> &grad, void *data);

  class SeldonOpt
  {

  private:

    nlopt_opt o;

    void mythrow(nlopt_result ret) const
    {
      switch (ret)
        {
        case NLOPT_FAILURE:
          throw Seldon::Error("SeldonOpt", "Nlopt failure.");
        case NLOPT_OUT_OF_MEMORY:
          throw Seldon::NoMemory("SeldonOpt",
                                 "Nlopt failed to allocate the"
                                 " requested storage space.");
        case NLOPT_INVALID_ARGS:
          throw Seldon::WrongArgument("SeldonOpt", "Nlopt invalid argument.");
        case NLOPT_ROUNDOFF_LIMITED:
          throw Seldon::Error("SeldonOpt", "Nlopt roundoff-limited.");
        case NLOPT_FORCED_STOP:
          throw Seldon::Error("SeldonOpt", "Nlopt forced stop.");
        default:
          break;
        }
    }


    typedef struct
    {
      SeldonOpt* o;
      mfunc mf;
      func f;
      void* f_data;
      svfunc vf;
      nlopt_munge munge_destroy, munge_copy;
    } myfunc_data;


    static void *free_myfunc_data(void *p)
    {
      myfunc_data *d = (myfunc_data *) p;
      if (d)
        {
          if (d->f_data && d->munge_destroy)
            d->munge_destroy(d->f_data);
          delete d;
        }
      return NULL;
    }


    static void *dup_myfunc_data(void *p)
    {
      myfunc_data *d = (myfunc_data *) p;
      if (d)
        {
          void *f_data;
          if (d->f_data && d->munge_copy) {
            f_data = d->munge_copy(d->f_data);
            if (!f_data)
              return NULL;
          }
          else
            f_data = d->f_data;
          myfunc_data *dnew = new myfunc_data;
          if (dnew)
            {
              *dnew = *d;
              dnew->f_data = f_data;
            }
          return (void*) dnew;
        }
      else
        return NULL;
    }


    static double myfunc(unsigned n, const double *x, double *grad,
                         void *d_)
    {
      myfunc_data *d = reinterpret_cast<myfunc_data*>(d_);
      return d->f(n, x, grad, d->f_data);
    }


    static void mymfunc(unsigned m, double *result,
                        unsigned n, const double *x,
                        double *grad, void *d_)
    {
      myfunc_data *d = reinterpret_cast<myfunc_data*>(d_);
      d->mf(m, result, n, x, grad, d->f_data);
      return;
    }


    Seldon::Vector<double> xtmp, gradtmp, gradtmp0;


    static double myvfunc(unsigned n, const double *x, double *grad,
                          void *d_)
    {
      myfunc_data *d = reinterpret_cast<myfunc_data*>(d_);
      Seldon::Vector<double> &xv = d->o->xtmp;
      if (n)
        memcpy(xv.GetData(), x, n * sizeof(double));

      double val=d->vf(xv, grad ? d->o->gradtmp : d->o->gradtmp0,
                       d->f_data);
      if (grad && n)
        {
          Seldon::Vector<double> &gradv = d->o->gradtmp;
          memcpy(grad, gradv.GetData(), n * sizeof(double));
        }
      return val;
    }


    void alloc_tmp()
    {
      if (xtmp.GetSize() != int(nlopt_get_dimension(o)))
        {
          xtmp = Seldon::Vector<double>(nlopt_get_dimension(o));
          gradtmp = Seldon::Vector<double>(nlopt_get_dimension(o));
        }
    }


    result last_result;
    double last_optf;
    nlopt_result forced_stop_reason;


  public:

    SeldonOpt() : o(NULL), xtmp(0), gradtmp(0), gradtmp0(0),
                  last_result(nlopt::FAILURE), last_optf(HUGE_VAL),
                  forced_stop_reason(NLOPT_FORCED_STOP)
    {
    }


    ~SeldonOpt()
    {
      nlopt_destroy(o);
    }


    SeldonOpt(algorithm a, unsigned n) :
      o(nlopt_create(nlopt_algorithm(a), n)),
      xtmp(0), gradtmp(0), gradtmp0(0),
      last_result(nlopt::FAILURE), last_optf(HUGE_VAL),
      forced_stop_reason(NLOPT_FORCED_STOP)
    {
      if (!o)
        throw Seldon::NoMemory("SeldonOpt(algorithm a, unsigned n)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      nlopt_set_munge(o, free_myfunc_data, dup_myfunc_data);
    }


    SeldonOpt(const SeldonOpt& f) :
      o(nlopt_copy(f.o)),
      xtmp(f.xtmp), gradtmp(f.gradtmp),
      gradtmp0(0),
      last_result(f.last_result),
      last_optf(f.last_optf),
      forced_stop_reason(f.forced_stop_reason)
    {
      if (f.o && !o)
        throw Seldon::NoMemory("SeldonOpt(const SeldonOpt& f)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
    }


    SeldonOpt& operator=(SeldonOpt const& f)
    {
      if (this == &f)
        return *this;
      nlopt_destroy(o);
      o = nlopt_copy(f.o);
      if (f.o && !o)
        throw Seldon::NoMemory("SeldonOpt::operator=(const SeldonOpt& f)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      xtmp = f.xtmp;
      gradtmp = f.gradtmp;
      last_result = f.last_result;
      last_optf = f.last_optf;
      forced_stop_reason = f.forced_stop_reason;
      return *this;
    }


    result optimize(Seldon::Vector<double>& x, double& opt_f)
    {
      if (o && int(nlopt_get_dimension(o)) != x.GetSize())
        throw Seldon::WrongArgument("SeldonOpt::optimize("
                                    "Seldon::Vector<double>& x,"
                                    " double& opt_f)", "Dimension mismatch.");
      forced_stop_reason = NLOPT_FORCED_STOP;
      nlopt_result ret =
        nlopt_optimize(o, x.GetSize() == 0 ? NULL : x.GetData(),
                       &opt_f);
      last_result = result(ret);
      last_optf = opt_f;
      if (ret == NLOPT_FORCED_STOP)
        mythrow(forced_stop_reason);
      mythrow(ret);
      return last_result;
    }


    Seldon::Vector<double> optimize(const Seldon::Vector<double> &x0)
    {
      Seldon::Vector<double> x(x0);
      last_result = optimize(x, last_optf);
      return x;
    }


    result last_optimize_result() const
    {
      return last_result;
    }


    double last_optimum_value() const
    {
      return last_optf;
    }


    algorithm get_algorithm() const
    {
      if (!o)
        throw Seldon::Error("SeldonOpt::get_algorithm() const",
                            "Uninitialized nlopt::SeldonOpt.");
      return algorithm(nlopt_get_algorithm(o));
    }


    const char *get_algorithm_name() const
    {
      if (!o)
        Seldon::Error("SeldonOpt::get_algorithm() const",
                      "Uninitialized nlopt::SeldonOpt.");
      return nlopt_algorithm_name(nlopt_get_algorithm(o));
    }


    unsigned get_dimension() const
    {
      if (!o)
        Seldon::Error("SeldonOpt::get_algorithm() const",
                      "Uninitialized nlopt::SeldonOpt.");
      return nlopt_get_dimension(o);
    }


    void set_min_objective(func f, void *f_data)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_min_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_set_min_objective(o, myfunc, d));
    }


    void set_min_objective(svfunc vf, void *f_data) {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_min_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = NULL;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = vf;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_set_min_objective(o, myvfunc, d));
      alloc_tmp();
    }


    void set_max_objective(func f, void *f_data) {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_min_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_set_max_objective(o, myfunc, d));
    }


    void set_max_objective(svfunc vf, void *f_data)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_max_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = NULL;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = vf;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_set_max_objective(o, myvfunc, d));
      alloc_tmp();
    }


    void set_min_objective(func f, void *f_data,
                           nlopt_munge md, nlopt_munge mc)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_min_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = md;
      d->munge_copy = mc;
      mythrow(nlopt_set_min_objective(o, myfunc, d));
    }


    void set_max_objective(func f, void *f_data,
                           nlopt_munge md, nlopt_munge mc)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::set_max_objective(func f, "
                               "void *f_data)", "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = md;
      d->munge_copy = mc;
      mythrow(nlopt_set_max_objective(o, myfunc, d));
    }


    void remove_inequality_constraints()
    {
      nlopt_result ret = nlopt_remove_inequality_constraints(o);
      mythrow(ret);
    }


    void add_inequality_constraint(func f, void *f_data, double tol = 0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_add_inequality_constraint(o, myfunc, d, tol));
    }


    void add_inequality_constraint(svfunc vf, void *f_data,
                                   double tol = 0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = NULL;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = vf;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_add_inequality_constraint(o, myvfunc, d, tol));
      alloc_tmp();
    }


    void add_inequality_mconstraint(mfunc mf, void *f_data,
                                    const Seldon::Vector<double> &tol)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_mconstraint(func f,"
                               " void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->mf = mf;
      d->f_data = f_data;
      d->f = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_add_inequality_mconstraint(o, tol.GetSize(),
                                               mymfunc, d,
                                               tol.GetSize() == 0
                                               ? NULL : tol.GetData()));
    }


    void remove_equality_constraints()
    {
      nlopt_result ret = nlopt_remove_equality_constraints(o);
      mythrow(ret);
    }


    void add_equality_constraint(func f, void *f_data, double tol = 0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_equality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      nlopt_add_equality_constraint(o, myfunc, d, tol);
      mythrow(nlopt_add_equality_constraint(o, myfunc, d, tol));
    }


    void add_equality_constraint(svfunc vf, void *f_data, double tol = 0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_equality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = NULL;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = vf;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_add_equality_constraint(o, myvfunc, d, tol));
      alloc_tmp();
    }


    void add_equality_mconstraint(mfunc mf, void *f_data,
                                  const Seldon::Vector<double> &tol)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_equality_mconstraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->mf = mf;
      d->f_data = f_data;
      d->f = NULL;
      d->vf = NULL;
      d->munge_destroy = d->munge_copy = NULL;
      mythrow(nlopt_add_equality_mconstraint(o, tol.GetSize(),
                                             mymfunc, d,
                                             tol.GetSize() == 0 ?
                                             NULL : tol.GetData()));
    }


    // For internal use in SWIG wrappers (see also above)
    void add_inequality_constraint(func f, void *f_data,
                                   nlopt_munge md, nlopt_munge mc,
                                   double tol=0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = md;
      d->munge_copy = mc;
      mythrow(nlopt_add_inequality_constraint(o, myfunc, d, tol));
    }


    void add_equality_constraint(func f, void *f_data,
                                 nlopt_munge md, nlopt_munge mc,
                                 double tol=0)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_constraint(func f, "
                               "void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->f = f;
      d->f_data = f_data;
      d->mf = NULL;
      d->vf = NULL;
      d->munge_destroy = md;
      d->munge_copy = mc;
      mythrow(nlopt_add_equality_constraint(o, myfunc, d, tol));
    }


    void add_inequality_mconstraint(mfunc mf, void *f_data,
                                    nlopt_munge md, nlopt_munge mc,
                                    const Seldon::Vector<double> &tol)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_inequality_mconstraint(func f,"
                               " void *f_data, double tol = 0)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->mf = mf;
      d->f_data = f_data;
      d->f = NULL;
      d->vf = NULL;
      d->munge_destroy = md; d->munge_copy = mc;
      mythrow(nlopt_add_inequality_mconstraint(o, tol.GetSize(),
                                               mymfunc, d,
                                               tol.GetSize() == 0
                                               ? NULL : tol.GetData()));
    }


    void add_equality_mconstraint(mfunc mf, void *f_data,
                                  nlopt_munge md, nlopt_munge mc,
                                  const Seldon::Vector<double> &tol)
    {
      myfunc_data *d = new myfunc_data;
      if (!d)
        throw Seldon::NoMemory("SeldonOpt::add_equality_mconstraint(mfunc mf,"
                               " void *f_data,"
                               "nlopt_munge md, nlopt_munge mc,"
                               "const Seldon::Vector<double> &tol)",
                               "Nlopt failed to allocate the"
                               " requested storage space.");
      d->o = this;
      d->mf = mf;
      d->f_data = f_data;
      d->f = NULL;
      d->vf = NULL;
      d->munge_destroy = md;
      d->munge_copy = mc;
      mythrow(nlopt_add_equality_mconstraint(o, tol.GetSize(), mymfunc,
                                             d, tol.GetSize() == 0
                                             ? NULL :
                                             tol.GetData()));
    }


#define SELDON_NLOPT_GETSET_VEC(name)                                   \
    void set_##name(double val) {                                       \
      mythrow(nlopt_set_##name##1(o, val));                             \
    }                                                                   \
    void get_##name(Seldon::Vector<double> &v) const {                  \
      if (o && int(nlopt_get_dimension(o)) != v.GetSize())              \
        throw Seldon::WrongArgument("SeldonOpt::get_" #name "(Vector&)" \
                                    " const",                           \
                                    "Nlopt invalid argument.");         \
      mythrow(nlopt_get_##name(o, v.GetSize() == 0 ? NULL :             \
                               v.GetData()));                           \
    }                                                                   \
    Seldon::Vector<double> get_##name() const {                         \
      if (!o) throw                                                     \
                Seldon::Error("SeldonOpt::get_" #name "() const",       \
                              "Uninitialized nlopt::SeldonOpt.");       \
      Seldon::Vector<double> v(nlopt_get_dimension(o));                 \
      get_##name(v);                                                    \
      return v;                                                         \
    }                                                                   \
    void set_##name(const Seldon::Vector<double> &v) {                  \
      if (o && int(nlopt_get_dimension(o)) != v.GetSize())              \
        throw Seldon::WrongArgument("SeldonOpt::get_" #name "(Vector&)" \
                                    " const",                           \
                                    "Nlopt invalid argument.");         \
      mythrow(nlopt_set_##name(o, v.GetSize() == 0 ?                    \
                               NULL : v.GetData()));                    \
    }

    SELDON_NLOPT_GETSET_VEC(lower_bounds)
    SELDON_NLOPT_GETSET_VEC(upper_bounds)


#define SELDON_NLOPT_GETSET(T, name)                                    \
      T get_##name() const {                                            \
        if (!o) throw                                                   \
                  Seldon::Error("SeldonOpt::get_" #name "() const",     \
                                "Uninitialized nlopt::SeldonOpt.");     \
        return nlopt_get_##name(o);                                     \
      }                                                                 \
      void set_##name(T name) {                                         \
        mythrow(nlopt_set_##name(o, name));                             \
      }


    SELDON_NLOPT_GETSET(double, stopval)
    SELDON_NLOPT_GETSET(double, ftol_rel)
    SELDON_NLOPT_GETSET(double, ftol_abs)
    SELDON_NLOPT_GETSET(double, xtol_rel)
    SELDON_NLOPT_GETSET_VEC(xtol_abs)
    SELDON_NLOPT_GETSET(int, maxeval)
    SELDON_NLOPT_GETSET(double, maxtime)

    SELDON_NLOPT_GETSET(int, force_stop)


    void force_stop()
    {
      set_force_stop(1);
    }


    void set_local_optimizer(const SeldonOpt &lo)
    {
      nlopt_result ret = nlopt_set_local_optimizer(o, lo.o);
      mythrow(ret);
    }


    SELDON_NLOPT_GETSET(unsigned, population)
    SELDON_NLOPT_GETSET_VEC(initial_step)


    void set_default_initial_step(const Seldon::Vector<double>& x)
    {
      nlopt_result ret
        = nlopt_set_default_initial_step(o, x.GetSize() == 0 ?
                                         NULL : x.GetData());
      mythrow(ret);
    }


    void get_initial_step(const Seldon::Vector<double> &x,
                          Seldon::Vector<double> &dx) const
    {
      if (o && (int(nlopt_get_dimension(o)) != x.GetSize()
                || int(nlopt_get_dimension(o)) != dx.GetSize()))
        throw Seldon::WrongArgument("SeldonOpt::get_initial_step("
                                    "Vector<double>& x, double& dx)",
                                    "Dimension mismatch.");
      nlopt_result ret = nlopt_get_initial_step(o, x.GetSize() == 0 ?
                                                NULL : x.GetData(),
                                                dx.GetSize() == 0 ?
                                                NULL : dx.GetData());
      mythrow(ret);
    }


    Seldon::Vector<double> get_initial_step_(const Seldon::Vector<double>& x)
    const
    {
      if (!o)
        throw Seldon::Error("SeldonOpt::get_initial_step_"
                            "(const Seldon::Vector<double>& x)",
                            "Uninitialized nlopt::SeldonOpt.");
      Seldon::Vector<double> v(nlopt_get_dimension(o));
      get_initial_step(x, v);
      return v;
    }
  };


#undef SELDON_NLOPT_GETSET
#undef SELDON_NLOPT_GETSET_VEC


}

#endif
