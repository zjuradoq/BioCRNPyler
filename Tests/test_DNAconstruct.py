
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler import *
import copy

def promoter_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    x = DNA_construct([P],mechanisms = mechs,parameters=parameters)
    assert(not(P in x.parts_list)) #make sure P is copied
    assert(P in x) #testing "contains" function of Construct
    y = x.enumerate_components()
    assert([x[0]] == y)
    assert(x[0].transcript is None)

    #works in reverse as well
    x = DNA_construct([[P,"reverse"]],mechanisms = mechs,parameters=parameters)
    assert(not(P in x.parts_list)) #make sure P is copied
    assert(P in x) #testing "contains" function of Construct
    y = x.enumerate_components()
    assert([x[0]] == y)
    assert(x[0].transcript is None)

def promoter_terminator_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    #minimal RNA transcription
    x = DNA_construct([P,T],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

    #"incoherent" transcription
    x = DNA_construct([[P,"reverse"],T],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y == [x[0]]) #correct promoter is returned

    #transcription in reverse works as well
    x = DNA_construct([T,[P,"reverse"]],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y[0] == x[1]) #correct promoter is returned
    assert(x[1].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[1]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

def circular_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    #circular construct
    x = DNA_construct([P,T],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs


    x = DNA_construct([[P,"reverse"],T],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

    x = DNA_construct([[P,"reverse"],[T,"reverse"],],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()

    z = DNA_construct([T,P],mechanisms = mechs,parameters=parameters,circular=True)
    e = z.enumerate_components()
    assert(y[1]==e[1]) #different DNA constructs lead to the same RNA construct
    assert(y[0]==x[0]) #correct promoter is working
    assert(e[0]==z[1]) #correct promoter is working