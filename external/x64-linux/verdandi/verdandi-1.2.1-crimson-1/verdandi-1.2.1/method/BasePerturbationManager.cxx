// Copyright (C) 2010-2011 INRIA
// Author(s): Anne Tilloy, Vivien Mallet, Kévin Charpentier
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_METHOD_BASEPERTURBATIONMANAGER_CXX


#include "BasePerturbationManager.hxx"
#include "seldon/vector/VectorCollection.cxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! Builds the manager. */
    template <class Derived>
    BasePerturbationManager<Derived>::BasePerturbationManager()
    {
        /*** Initializations ***/

        MessageHandler::AddRecipient("perturbation_manager", *this,
                                     BasePerturbationManager::StaticMessage);
    }


    //! Main constructor.
    /*! Builds the manager.
      \param[in] configuration_file configuration file.
    */
    template <class Derived>
    BasePerturbationManager<Derived>
    ::BasePerturbationManager(string configuration_file)
    {
        MessageHandler::AddRecipient("perturbation_manager", *this,
                                     BasePerturbationManager::StaticMessage);
    }


    //! Destructor.
    template <class Derived>
    BasePerturbationManager<Derived>::~BasePerturbationManager()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the manager.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class Derived>
    void BasePerturbationManager<Derived>
    ::Initialize(string configuration_file)
    {
    }


    /*! \brief Generates a random vector or several vectors according to some
      distribution. */
    /*! This method generates one vector or several vectors, according to the
      shape of \a output on entry. If \a variance is a \f$ N \times N \f$
      matrix, then \a output must be a vector of size \f$ m N \f$, and this
      method generates a sample of \f$ m \f$ vectors of size \f$ N \f$, stored
      in \a output. The sample can be generated with or without correlation
      between the vectors, depending on \a correlation.

      \param[in] pdf probability density function: "Normal" or "LogNormal".
      \param[in] variance covariance matrix of the distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal (log-normal) distribution, any
      component \f$ i \f$ of the vector (the logarithm of the vector) lies in
      \f$ [\mu_i - a \sigma_i, \mu_i + b \sigma_i] \f$ where \f$ \mu_i \f$ is
      the mean (median) of the random component \f$ i \f$ and \f$ \sigma_i \f$
      is its standard deviation (the standard deviation of its logarithm).
      \param[in] correlation if non-empty, correlations with the first
      vector: \a correlation should then contain \f$ m - 1 \f$ elements if
      there are \f$ m \f$ vectors to be generated.
      \param[in,out] output on entry, mean or median vector(s); on exit, the
      sample. If the size of \a output is not a multiple of \f$ N \f$, an
      exception is thrown.
    */
    template<class Derived>
    template <class T0, class Prop0, class Allocator0,
              class T1, class Allocator1>
    void BasePerturbationManager<Derived>
    ::Sample(string pdf,
             Matrix<T0, Prop0, RowSymPacked, Allocator0>& variance,
             Vector<double, VectFull>& parameter,
             Vector<double, VectFull>& correlation,
             Vector<T1, VectFull, Allocator1>& output)
    {
        if (variance.GetM() != variance.GetN())
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<VectFull>)",
                                "The covariance matrix must be a square "
                                "matrix, but a " + to_str(variance.GetM())
                                + " x " + to_str(variance.GetN())
                                + " matrix was provided.");

        int N = variance.GetN();
        int Nvector = output.GetLength() / N;

        if (output.GetLength() % N != 0)
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<VectFull>)",
                                "The " + to_str(variance.GetM())
                                + " x " + to_str(variance.GetN())
                                + " covariance matrix has dimensions "
                                + "incompatible with the mean(s) or "
                                + "median(s) whose cumulated size is "
                                + to_str(output.GetLength()) + ".");

        if (correlation.GetLength() != 0
            && Nvector - 1 != correlation.GetLength())
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<VectFull>)",
                                to_str(Nvector) + " vector(s) is/are "
                                + "provided, but "
                                + to_str(correlation.GetLength())
                                + " correlation(s) is/are given.");

        if(pdf != "Normal" && pdf != "LogNormal")
            throw ErrorArgument("BasePerturbationManager"
                                "::Sample(string, Matrix, Vector,"
                                " Vector, Vector)",
                                "The distribution \""
                                + pdf + "\" is not supported.");

        if (Nvector == 1 && pdf == "Normal")
            static_cast<Derived*>(this)->Normal(variance,
                                                parameter, output);

        else if (pdf == "Normal")
        {
            Vector<T1, VectFull, Allocator1> first_vector;
            first_vector.SetData(N, &output(0));

            static_cast<Derived*>(this)->Normal(variance,
                                                parameter, first_vector);

            for (int i = 1; i < Nvector; i++)
            {
                Vector<T1, VectFull, Allocator1> vector_i;
                vector_i.SetData(N, &output(i * N));

                if (correlation.GetLength() != 0 && correlation(i - 1) == 1.)
                    Copy(first_vector, vector_i);
                else
                {
                    static_cast<Derived*>(this)->Normal(variance, parameter,
                                                        vector_i);
                    if (correlation.GetLength() != 0)
                    {
                        Mlt(T1(sqrt(1. - correlation(i-1)
                                    * correlation(i-1))), vector_i);
                        Add(T1(correlation(i - 1)), first_vector, vector_i);
                    }
                }
                vector_i.Nullify();
            }
            first_vector.Nullify();
        }

        else if (Nvector == 1 && pdf == "LogNormal")
            static_cast<Derived*>(this)->LogNormal(variance,
                                                   parameter, output);

        else if (pdf == "LogNormal")
        {
            Vector<T1, VectFull, Allocator1> perturbation(output.GetLength());
            perturbation.Fill(T1(0));

            Vector<T1, VectFull, Allocator1> first_vector;
            first_vector.SetData(N, &perturbation(0));

            static_cast<Derived*>(this)->Normal(variance,
                                                parameter, first_vector);

            for (int i = 1; i < Nvector; i++)
            {
                Vector<T1, VectFull, Allocator1> vector_i;
                vector_i.SetData(N, &perturbation(i * N));

                if (correlation.GetLength() != 0 && correlation(i - 1) == 1.)
                    Copy(first_vector, vector_i);
                else
                {
                    static_cast<Derived*>(this)->Normal(variance,
                                                        parameter,
                                                        vector_i);

                    if (correlation.GetLength() != 0)
                    {
                        Mlt(T1(sqrt(1. - correlation(i-1)*correlation(i-1))),
                            vector_i);
                        Add(T1(correlation(i - 1)), first_vector, vector_i);
                    }
                }
                vector_i.Nullify();
            }
            first_vector.Nullify();

            for (int k = 0; k < perturbation.GetM(); k++)
                output(k) *= exp(perturbation(k));
        }
    }


    //! Generates a random vector collection according to some distribution.
    /*!
      \param[in] pdf probability density function: "Normal" or "LogNormal".
      \param[in] variance covariance matrix of the distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal (log-normal) distribution, any
      component \f$ i \f$ of the vector (the logarithm of the vector) lies in
      \f$ [\mu_i - a \sigma_i, \mu_i + b \sigma_i] \f$ where \f$ \mu_i \f$ is
      the mean (median) of the random component \f$ i \f$ and \f$ \sigma_i \f$
      is its standard deviation (the standard deviation of its logarithm).
      \param[in] correlation if non-empty, correlations with the first
      vector: \a correlation should then contain \f$ m - 1 \f$ elements if
      there are \f$ m \f$ vectors to be generated.
      \param[in,out] output on entry, mean or median vectors; on exit, the
      sample in the form of a vector collection.
    */
    template <class Derived>
    template <class T0, class Prop0, class Allocator0,
              class T1, class Allocator1>
    void BasePerturbationManager<Derived>
    ::Sample(string pdf,
             Matrix<T0, Prop0, RowSymPacked, Allocator0>& variance,
             Vector<double, VectFull>& parameter,
             Vector<double, VectFull>& correlation,
             Vector<T1, Collection, Allocator1>& output)
    {

        if (variance.GetM() != variance.GetN())
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<Collection>)",
                                "The covariance matrix must be a square "
                                "matrix, but a " + to_str(variance.GetM())
                                + " x " + to_str(variance.GetN())
                                + " matrix was provided.");

        int N = variance.GetN();
        int Nvector = output.GetNvector();

        for (int i = 0; i < Nvector; i++)
            if (variance.GetN() != output.GetVector(i).GetLength())
                throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                    "Vector<Collection>)",
                                    "The " + to_str(variance.GetM())
                                    + " x " + to_str(variance.GetN())
                                    + " covariance matrix has dimensions "
                                    + "incompatible with the " + to_str(i)
                                    + "-th mean or median which size is "
                                    + to_str(output.GetVector(i).GetLength())
                                    + ".");

        if (correlation.GetLength() != 0
            && Nvector - 1 != correlation.GetLength())
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<Collection>)",
                                to_str(Nvector) + " vector(s) is/are "
                                + "provided, but "
                                + to_str(correlation.GetLength())
                                + " correlation(s) is/are given.");


        if(pdf != "Normal" && pdf != "LogNormal")
            throw ErrorArgument("BasePerturbationManager"
                                "::Sample(string, Matrix, Vector,"
                                " Vector, VectorCollection)",
                                "The distribution \""
                                + pdf + "\" is not supported.");

        if (Nvector == 1 && pdf == "Normal")
            static_cast<Derived*>(this)-> Normal(variance, parameter,
                                                 output.GetVector(0));

        else if (pdf == "Normal")
        {

            static_cast<Derived*>(this)->Normal(variance, parameter,
                                                output.GetVector(0));

            for (int i = 1; i < Nvector; i++)
            {
                if (correlation.GetLength() != 0 && correlation(i - 1) == 1.)
                    output.GetVector(i) = output.GetVector(0);
                else
                {
                    static_cast<Derived*>(this)->
                        Normal(variance,
                               parameter,
                               output.GetVector(i));

                    if (correlation.GetLength() != 0)
                    {
                        Mlt(typename T1::value_type(
                                sqrt(1. -
                                     correlation(i-1) * correlation(i-1))),
                            output.GetVector(i));
                        Add(typename T1::value_type(correlation(i - 1)),
                            output.GetVector(0),
                            output.GetVector(i));
                    }
                }
            }
        }

        else if (Nvector == 1 && pdf == "LogNormal")
            static_cast<Derived*>(this)->
                LogNormal(variance, parameter,
                          output.GetVector(0));


        else if (pdf == "LogNormal")
        {
            Vector<T1, Collection, Allocator1> perturbation(Nvector);
            T1* tmp;
            for (int i = 0; i < Nvector; i++)
            {
                tmp = new T1(N);
                tmp->Fill(typename T1::value_type(0));
                perturbation.SetVector(i, *tmp);
                tmp->Nullify();
                delete tmp;
            }

            static_cast<Derived*>(this)->
                Normal(variance,
                       parameter,
                       perturbation.GetVector(0));

            for (int i = 1; i < Nvector; i++)
            {
                if (correlation.GetLength() != 0 && correlation(i - 1) == 1.)
                    perturbation.GetVector(i) = perturbation.GetVector(0);
                else
                {
                    static_cast<Derived*>(this)->
                        Normal(variance, parameter,
                               perturbation.GetVector(i));

                    if (correlation.GetLength() != 0)
                    {
                        Mlt(typename T1::value_type(
                                sqrt(1. -
                                     correlation(i-1) * correlation(i-1))),
                            perturbation.GetVector(i));
                        Add(typename T1::value_type(correlation(i - 1)),
                            perturbation.GetVector(0),
                            perturbation.GetVector(i));
                    }
                }
            }

            for (int k = 0; k < perturbation.GetM(); k++)
                output(k) *= exp(perturbation(k));
            perturbation.Deallocate();
        }
    }


    /*! \brief Generates a random vector or several vectors according to a
      homogeneous distribution. */
    /*! This method generates one vector or several vectors, according to the
      shape of \a output on entry. If \a correlation is of size \f$ m - 1 \f$,
      then \a output must be a vector of size \f$ m N \f$, and this method
      generates a sample of \f$ m \f$ vectors of size \f$ N \f$, stored in \a
      output. The sample can be generated with or without correlation between
      the vectors, depending on \a correlation.

      For a normal (log-normal) homogeneous distribution, for each vector of
      size \f$ N \f$, the same random number, with centered normal
      distribution and variance \a variance, is added to every component of
      the vector (logarithm of the vector).

      \param[in] pdf probability density function: "NormalHomogeneous" or
      "LogNormalHomogeneous".
      \param[in] variance variance of the distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal (log-normal) distribution, any random
      value lies in \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$
      is the mean (median) of the random variable and \f$ \sigma \f$ is its
      standard deviation (the standard deviation of its logarithm).
      \param[in] correlation if non-empty, correlations with the first
      vector: \a correlation should then contain \f$ m - 1 \f$ elements if
      there are \f$ m \f$ vectors to be sampled.
      \param[in,out] output on entry, mean or median vector(s); on exit, the
      sample. If the size of \a output is not a multiple of that of \a
      correlation plus one, an exception is thrown.
    */
    template <class Derived>
    template <class T0,
              class T1, class Allocator1>
    void BasePerturbationManager<Derived>
    ::Sample(string pdf, T0 variance,
             Vector<double, VectFull>& parameter,
             Vector<double, VectFull>& correlation,
             Vector<T1, VectFull, Allocator1>& output)
    {
        if (output.GetLength() % (correlation.GetLength() + 1) != 0)
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<VectFull>)",
                                "The size of the mean(s) or median(s) ("
                                + to_str(output.GetLength()) + ") is "
                                + "incompatible with the size of the "
                                + "correlation vector ("
                                + to_str(correlation.GetLength()) + ").");

        int Nvector = correlation.GetLength() + 1;
        int N = output.GetLength() / Nvector;

        if(pdf != "NormalHomogeneous" && pdf != "LogNormalHomogeneous")
            throw ErrorArgument("BasePerturbationManager"
                                "::Sample(string, double, Vector,"
                                " Vector, VectorCollection)",
                                "The distribution \""
                                + pdf + "\" is not supported.");

        if (Nvector == 1 && pdf == "NormalHomogeneous")
            static_cast<Derived*>(this)->NormalHomogeneous(variance,
                                                           parameter, output);

        else if (pdf == "NormalHomogeneous")
        {

            Vector<T1, VectFull, Allocator1> first_vector(N), vector_i(N);
            first_vector.SetData(N, &output(0));

            static_cast<Derived*>(this)->NormalHomogeneous(variance,
                                                           parameter,
                                                           first_vector);

            for (int i = 1; i < Nvector; i++)
            {
                vector_i.SetData(N, &output(i * N));
                static_cast<Derived*>(this)->NormalHomogeneous(variance,
                                                               parameter,
                                                               vector_i);
                Mlt(T1(1. - correlation(i - 1)), vector_i);
                Add(T1(correlation(i - 1)), first_vector, vector_i);
                vector_i.Nullify();
            }
            first_vector.Nullify();
        }

        else if (Nvector == 1 && pdf == "LogNormalHomogeneous")
            static_cast<Derived*>(this)->LogNormalHomogeneous(variance,
                                                              parameter,
                                                              output);

        else if (pdf == "LogNormalHomogeneous")
        {
            Vector<T1, VectFull, Allocator1> perturbation(output.GetLength());
            perturbation.Fill(T1(0));

            Vector<T1, VectFull, Allocator1> first_vector(N), vector_i(N);
            first_vector.SetData(N, &perturbation(0));

            static_cast<Derived*>(this)->NormalHomogeneous(variance,
                                                           parameter,
                                                           first_vector);

            for (int i = 1; i < Nvector; i++)
            {
                vector_i.SetData(N, &perturbation(i * N));
                static_cast<Derived*>(this)->NormalHomogeneous(variance,
                                                               parameter,
                                                               vector_i);
                Mlt(T1(sqrt(1. - correlation(i - 1) * correlation(i - 1))),
                    vector_i);
                Add(T1(correlation(i - 1)), first_vector, vector_i);
                vector_i.Nullify();
            }
            first_vector.Nullify();

            for (int k = 0; k < perturbation.GetLength(); k++)
                output(k) *= exp(perturbation(k));
        }
    }


    /*! \brief Generates a random vector collection according to a homogeneous
      distribution. */
    /*! The number of vectors in \a output should be the size of \a
      correlation plus one.

      For a normal (log-normal) homogeneous distribution, for each vector, the
      same random number, with centered normal distribution and variance \a
      variance, is added to every component of the vector (logarithm of the
      vector).

      \param[in] pdf probability density function: "NormalHomogeneous" or
      "LogNormalHomogeneous".
      \param[in] variance variance of the distribution.
      \param[in] parameter vector of parameters. The vector may either be
      empty or contain two clipping parameters \f$ (a, b) \f$. With the
      clipping parameters, for a normal (log-normal) distribution, any random
      value lies in \f$ [\mu - a \sigma, \mu + b \sigma] \f$ where \f$ \mu \f$
      is the mean (median) of the random variable and \f$ \sigma \f$ is its
      standard deviation (the standard deviation of its logarithm).
      \param[in] correlation if non-empty, correlations with the first
      vector: \a correlation should then contain \f$ m - 1 \f$ elements if
      there are \f$ m \f$ vectors to be generated.
      \param[in,out] output output on entry, mean or median vectors of the
      vector collection; on exit, the sample in the form of a vector
      collection.
    */
    template <class Derived>
    template <class T0,
              class T1, class Allocator1>
    void BasePerturbationManager<Derived>
    ::Sample(string pdf, T0 variance,
             Vector<double, VectFull>& parameter,
             Vector<double, VectFull>& correlation,
             Vector<T1, Collection, Allocator1>& output)
    {
        int Nvector = output.GetNvector();
        int N = output.GetVector(0).GetLength();

        for (int i = 0; i < Nvector; i++)
            if (output.GetVector(i).GetLength() != N)
                throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                    "Vector<Collection>)",
                                    + "The size of the " + to_str(i) +
                                    + "-th mean or median ("
                                    + to_str(output.GetVector(i).GetLength())
                                    + ") is incompatible with that of the "
                                    + "first mean or median ("
                                    + to_str(N) + ").");

        if (correlation.GetLength() != 0
            && Nvector - 1 != correlation.GetLength())
            throw ErrorArgument("BasePerturbationManager::Sample(..., "
                                "Vector<Collection>)",
                                to_str(Nvector) + " vector(s) is/are "
                                + "provided, but "
                                + to_str(correlation.GetLength())
                                + " correlation(s) is/are given.");

        if(pdf != "NormalHomogeneous" && pdf != "LogNormalHomogeneous")
            throw ErrorArgument("BasePerturbationManager"
                                "::Sample(string, double, Vector,"
                                " Vector, VectorCollection)",
                                "The distribution \""
                                + pdf + "\" is not supported.");

        if (Nvector == 1 && pdf == "NormalHomogeneous")
            static_cast<Derived*>(this)->
                NormalHomogeneous(variance,
                                  parameter,
                                  output.GetVector(0));

        else if (pdf == "NormalHomogeneous")
        {
            static_cast<Derived*>(this)->
                NormalHomogeneous(variance,
                                  parameter,
                                  output.GetVector(0));

            for (int i = 1; i < Nvector; i++)
            {
                static_cast<Derived*>(this)->
                    NormalHomogeneous(variance,
                                      parameter,
                                      output.GetVector(i));

                if (correlation.GetLength() != 0)
                {
                    Mlt(typename T1::value_type(
                            sqrt(1. - correlation(i-1) * correlation(i-1))),
                        output.GetVector(i));
                    Add(typename T1::value_type(correlation(i - 1)),
                        output.GetVector(0),
                        output.GetVector(i));
                }
            }
        }

        else if (Nvector == 1 && pdf == "LogNormalHomogeneous")
            static_cast<Derived*>(this)->
                LogNormalHomogeneous(variance,
                                     parameter,
                                     output.GetVector(0));

        else if (pdf == "LogNormalHomogeneous")
        {
            Vector<T1, Collection, Allocator1> perturbation(Nvector);
            T1* tmp;
            for (int i = 0; i < Nvector; i++)
            {
                tmp = new T1(N);
                tmp->Fill(typename T1::value_type(0));
                perturbation.SetVector(i, *tmp);
                tmp->Nullify();
                delete tmp;
            }

            static_cast<Derived*>(this)->
                NormalHomogeneous(variance,
                                  parameter,
                                  perturbation.GetVector(0));

            for (int i = 1; i < Nvector; i++)
            {
                static_cast<Derived*>(this)->
                    NormalHomogeneous(variance,
                                      parameter,
                                      perturbation.GetVector(i));

                if (correlation.GetLength() != 0)
                {
                    Mlt(sqrt(1. - correlation(i - 1) * correlation(i - 1)),
                        perturbation.GetVector(i));
                    Add(correlation(i - 1), perturbation.GetVector(0),
                        perturbation.GetVector(i));
                }
            }

            for (int k = 0; k < perturbation.GetM(); k++)
                output(k) *= exp(perturbation(k));
            perturbation.Deallocate();
        }
    }

} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BASEPERTURBATIONMANAGER_CXX
#endif
