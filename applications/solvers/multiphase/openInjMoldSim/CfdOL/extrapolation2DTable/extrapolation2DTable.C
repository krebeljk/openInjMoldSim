/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "openFoamTableReader.H"
#include "Vector.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::extrapolation2DTable<Type>::readTable()
{
    fileName fName(fileName_);
    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorIn
        (
            "Foam::extrapolation2DTable<Type>::readTable()"
        )   << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are in ascending order
    checkOrder();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable()
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(extrapolation2DTable::WARN),
    fileName_("fileNameIsUndefined"),
    reader_(NULL)
{}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable
(
    const List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >& values,
    const boundsHandling bounds,
    const fileName& fName
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(values),
    boundsHandling_(bounds),
    fileName_(fName),
    reader_(NULL)
{}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable(const fileName& fName)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(extrapolation2DTable::WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>(dictionary()))
{
    readTable();
}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable(const dictionary& dict)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(),
    boundsHandling_(wordToBoundsHandling(dict.lookup("outOfBounds"))),
    fileName_(dict.lookup("fileName")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
Foam::extrapolation2DTable<Type>::extrapolation2DTable
(
     const extrapolation2DTable& extrapTable
)
:
    List<Tuple2<scalar, List<Tuple2<scalar, Type> > > >(extrapTable),
    boundsHandling_(extrapTable.boundsHandling_),
    fileName_(extrapTable.fileName_),
    reader_(extrapTable.reader_)    // note: steals reader. Used in write().
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::extrapolation2DTable<Type>::extrapolateValue
(
    const List<Tuple2<scalar, Type> >& data,
    const scalar lookupValue
) const
{
    label n = data.size();

    scalar minLimit = data.first().first();
    scalar maxLimit = data.last().first();
    if (n == 1)
      return data.first().second();
    if (lookupValue < minLimit)
    {
        switch (boundsHandling_)
        {
            case extrapolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case extrapolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")" << nl
                    << "    extrapolating the first entry"
                    << endl;
                // fall-through to 'EXTRAPOLATE'
            }
            case extrapolation2DTable::EXTRAPOLATE:
            {
		scalar x1 = minLimit;
		scalar x2 = data[1].first();
		Type y1 = data[0].second();
		Type y2 = data[1].second();
                return y1 + (lookupValue - x1)/(x2-x1) * (y2 - y1); //extrapolation
                break;
            }
        }
    }
    else if (lookupValue > maxLimit)
    {
        switch (boundsHandling_)
        {
            case extrapolation2DTable::ERROR:
            {
                FatalErrorIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << exit(FatalError);
                break;
            }
            case extrapolation2DTable::WARN:
            {
                WarningIn
                (
                    "Foam::extrapolation2DTable<Type>::extrapolateValue"
                    "("
                        "List<Tuple2<scalar, Type> >&, "
                        "const scalar"
                    ")"
                )   << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")" << nl
                    << "    Continuing with the last entry"
                    << endl;
                // fall-through to 'EXTRAPOLATE'
            }
            case extrapolation2DTable::EXTRAPOLATE:
            {
		 scalar x1 = maxLimit;
		 scalar x2 = data[data.size()-2].first();
		 Type y1 = data[data.size()-1].second();
		 Type y2 = data[data.size()-2].second();
                return y1 + (lookupValue - x1)/(x2-x1) * (y2 - y1); //extrapolation
                return data.last().second();
                break;
            }
        }
    }

    // look for the correct range in X
    label lo = 0;
    label hi = 0;

    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= data[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return data[lo].second();
    }
    else
    {
        Type m =
            (data[hi].second() - data[lo].second())
           /(data[hi].first() - data[lo].first());

        // normal interpolation
        return data[lo].second() + m*(lookupValue - data[lo].first());
    }
}


template<class Type>
template<class L>
Foam::label Foam::extrapolation2DTable<Type>::Xi
(  
    const List<Tuple2<scalar, L> >& t,
    const scalar valueX,
    const bool reverse
) const
{
    label nX = t.size();
    //Binary or interpolation search for better efficience
    //First, cover the cases where the point lies outwards
    if(reverse)
    {
	if(t[nX - 1].first() <= valueX)
	  return nX - 1; //point lies outwards, return last index
	else if(t[0].first() > valueX)
	  return 1; //point lies outwards, return 2nd index
    }
    else
    {
	if(t[nX - 1].first() <= valueX)
	  return nX - 2; //point lies outwards, return second last index
	else if(t[0].first() > valueX)
	  return 0; //point lies outwards, return first index
    }
    
    label xmin = 0;
    label xmax = nX - 1;
    //label xtest = (xmax + xmin) / 2; //binary search
    label xtest = min(xmax - 1, max(xmin + 1, xmin + (valueX - t[xmin].first()) / (t[xmax].first() - t[xmin].first()) * (xmax - xmin))); //interpolation search
    while(xmax - xmin > 1)
    {
      if(t[xtest].first() <= valueX)
	xmin = xtest;
      else if(t[xtest].first() > valueX)
	xmax = xtest;
      if(t[xtest].first() == valueX)
	break;
      //xtest = (xmax + xmin) / 2; //binary search
      xtest = min(xmax - 1, max(xmin + 1, xmin + (valueX - t[xmin].first()) / (t[xmax].first() - t[xmin].first()) * (xmax - xmin))); //interpolation search
    }
    if(reverse)
      return xtest + 1;
    else
      return xtest;
    
    /*
    // Linear search. Probably not so efficient :D
    if (reverse)
    {
	for(int i = nX - 1; i > 0; i--)
	  if(t[i].first() <= valueX)
	  {
	    if(i == nX - 1) //point lies outwards, return last index
	      return nX - 1;
	    else
	      return i + 1;
	  }
	return 1; //point lies outwards, return 2nd index
    }
    else
    {
	for(int i = 0; i < nX - 1; i++)
	  if(t[i].first() > valueX)
	  {
	    if(i == 0)
	      return 0; //point lies outwards, return first index
	    else
	      return i-1;
	  }
	return nX - 2; //point lies outwards, return second last index
    }
    */
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Type Foam::extrapolation2DTable<Type>::operator()
(
    const scalar valueX,
    const scalar valueY
) const
{
    
    // Considers all of the list in Y being equal
    label nX = this->size();
    label nY = this->first().second().size();
    const table& t = *this;
    if (nX == 0)
    {
        WarningIn
        (
            "Type Foam::extrapolation2DTable<Type>::operator()"
            "("
                "const scalar, "
                "const scalar"
            ") const"
        )
            << "cannot extrapolate a zero-sized table - returning zero" << endl;

        return pTraits<Type>::zero;
    }
    else if (nX == 1)
    {
        // only 1 row (in X) - extrapolate to find Y value
        return extrapolateValue(t.first().second(), valueY);
    }
    else
    {
        // have 2-D data, extrapolate
	// find low and high indices in the X range that bound valueX
	label x0i = Xi(t, valueX, false);
	label x1i = x0i + 1; //Xi(t, valueX, true);
	const List<Tuple2<scalar, Type> >& row0 = t[x0i].second();
	const List<Tuple2<scalar, Type> >& row1 = t[x1i].second();
	//factor for interpolating between both rows
	scalar factor = (valueX - t[x0i].first()) / (t[x1i].first() - t[x0i].first());
	if(nY == 1)
	  return row0.first().second() + factor * (row1.first().second() - row0.first().second());
	// find low and high indices in the Y range that bound valueY
	label y00 = Xi(row0, valueY, false);
	label y01 = y00 + 1; //Xi(row0, valueY, true);
	label y10 = Xi(row1, valueY, false);
	label y11 = y10 + 1; //Xi(row1, valueY, true);
	
	label ymin = min(y00, y10);
	label ymax = max(y01, y11);
	List<vector> points;
	label index = 0;
	for(label i = ymin; i <= ymax; i++)
	{
	  points.append(vector(
	    valueX,
	    row0[i].first()  + factor * (row1[i].first()  - row0[i].first() ),
	    row0[i].second() + factor * (row1[i].second() - row0[i].second())
	  ));
	  if(points[i - ymin].y() < valueY)
	    index = i - ymin;
	}
	return (points[index] + (valueY - points[index].y()) / (points[index + 1].y() - points[index].y()) * (points[index + 1] - points[index])).z();
    }
}


template<class Type>
Foam::word Foam::extrapolation2DTable<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case extrapolation2DTable::ERROR:
        {
            enumName = "error";
            break;
        }
        case extrapolation2DTable::WARN:
        {
            enumName = "warn";
            break;
        }
        case extrapolation2DTable::EXTRAPOLATE:
        {
            enumName = "extrapolate";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::extrapolation2DTable<Type>::boundsHandling
Foam::extrapolation2DTable<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return extrapolation2DTable::ERROR;
    }
    else if (bound == "warn")
    {
        return extrapolation2DTable::WARN;
    }
    else if (bound == "extrapolate")
    {
        return extrapolation2DTable::EXTRAPOLATE;
    }
    else
    {
        WarningIn
        (
            "Foam::extrapolation2DTable<Type>::wordToBoundsHandling"
            "("
            "    const word&"
            ")"
        )   << "bad outOfBounds specifier " << bound << " using 'warn'" << endl;

        return extrapolation2DTable::WARN;
    }
}


template<class Type>
typename Foam::extrapolation2DTable<Type>::boundsHandling
Foam::extrapolation2DTable<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;
    return prev;
}


template<class Type>
void Foam::extrapolation2DTable<Type>::checkOrder() const
{
    label n = this->size();
    const table& t = *this;

    scalar prevValue = t[0].first();

    for (label i=1; i<n; ++i)
    {
        const scalar currValue = t[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorIn
            (
                "Foam::extrapolation2DTable<Type>::checkOrder() const"
            )   << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::extrapolation2DTable<Type>::write(Ostream& os) const
{
    os.writeKeyword("fileName")
        << fileName_ << token::END_STATEMENT << nl;
    os.writeKeyword("outOfBounds")
        << boundsHandlingToWord(boundsHandling_) << token::END_STATEMENT << nl;

    *this >> os;
}


// ************************************************************************* //
