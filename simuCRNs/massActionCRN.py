import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('dark_background')
from functools import reduce
import re
from collections import OrderedDict
from ipywidgets import interact, interactive, fixed, interact_manual, interactive_output
import scipy.integrate as integrate
import ipywidgets as widgets
import copy


class massActionJSON( object ):
    """ Parser for a mass action CRN specified in JSON.
        It will construct the object necessary to create a `massActionCRN` instanceself.

        Example for the following CRN spec:

        {
            "type": "mass action",

            "species": {
                "X": "0.5",
                "Y": "0.25",
                "Z": "0.25"
            },

            "reactions": {
                "2X + 3Y <-> Z": [ "0.0001", "0.000008" ],
                "2X -> Y + Z": [ "0.003" ],
                "X -> ": [ "0.229" ],
                "Z  + 5X <-> Y": [ "0.4", "0.001" ]
            }
        }

        It will create:

        species_name: ['X', 'Y', 'Z']
        x0: OrderedDict([('X', 0.5), ('Y', 0.25), ('Z', 0.25)])
        rates: OrderedDict([('k1', 0.0001), ('k-1', 8e-06), ('k2', 0.003), ('k3', 0.229), ('k4', 0.4), ('k-4', 0.001)])
        reactions:
        2X + 3Y <-> Z
        	 (array([-2., -3.,  0.]), array([0., 0., 1.]))
        	 (array([-0., -0., -1.]), array([ 2.,  3., -0.]))
        2X -> Y + Z
        	 (array([-2.,  0.,  0.]), array([0., 1., 1.]))
        X ->
        	 (array([-1.,  0.,  0.]), array([0., 0., 0.]))
        Z  + 5X <-> Y
        	 (array([-5.,  0., -1.]), array([0., 1., 0.]))
        	 (array([-0., -1., -0.]), array([ 5., -0.,  1.]))

    """
    def __init__( self, json_dict, debug = 0 ):
        self._json_dict = json_dict
        self._debug = debug
        self._debug_name = "[JSON MA PARSER]"

    def parse( self ):
        """ Parses the mass action JSON specification.
            Returns the corresponding CRN object.
        """
        self._name = self._json_dict['name']

        self._species_name = list( self._json_dict[ 'species' ].keys() )
        self._x0 = OrderedDict(map(lambda kv: (kv[0], float(kv[1])),
                               self._json_dict[ 'species' ].items()))

        self._rates = self._get_rates()
        self._reactions = self._get_reactions()

        if self._debug > 0:
            print( self._debug_name+'species:', self._species_name )
            print( self._debug_name+'x0:', self._x0 )
            print( self._debug_name+'rates:', self._rates )
            print( self._debug_name+'reactions:')
            for r in self._reactions:
                print( r )
                for rr in self._reactions[ r ]:
                    print("\t", rr)

        return massActionCRN( self._name, self._species_name, self._x0,
                              self._rates,self._reactions, self._debug-1 )

    def _get_rates( self ):
        """ Returns the rates in the appropriate format (c.f. class massActionCRN).
        """
        rates = OrderedDict( {} )
        for i, reaction in enumerate( self._json_dict[ 'reactions' ] ):
            rate_list = self._json_dict[ 'reactions'][ reaction ]
            for j, rate in enumerate( rate_list ):
                modif = "" if j == 0 else "-"
                rate_name = "k{}{}".format( modif, i+1 )
                rates[ rate_name ] = float( rate )
        return rates

    def _parse_reaction( self, reaction_str ):
        """ Transforms a reaction string "X + Y -> Z" in the pair
            (np.array([-1,-1,0]), np.array([0,0,1])).

            Returns a list of such pairs containing one or two pairs
            depending on if the reaction is written with "->" or "<->".
        """
        reaction_list = []

        clean_reaction = reaction_str.replace(" ", "")

        if '<->' in reaction_str:
            separator = '<->'
        else:
            separator = '->'

        central_split = clean_reaction.split( separator )

        reaction = []
        for i_side, side in enumerate( central_split ):
            species_split = side.split("+")

            reaction_side = np.zeros( len( self._x0 ) )

            if side != '':
                for s in species_split:
                    name, coeff = self._get_specie_and_stoc( s )
                    index_in_reaction_side = self._species_name.index( name )
                    reaction_side[ index_in_reaction_side ] = -1*coeff if i_side == 0 else coeff

            reaction.append( reaction_side )

        reaction_list.append( tuple( reaction ) )

        if separator == '<->':
            reaction_list.append( ( -1*reaction[ 1 ], -1*reaction[ 0 ] ) )

        return reaction_list

    def _get_specie_and_stoc( self, specie_in_reac ):
        """
            specie_in_reac: string of the form <float><specie_name> e.g.: 2X, 2.89Y_0

            Warning: the function will raise exception if:
                - no specie is found in the string e.g.: '', '234'
                - <specie_name> was not listed in the list of species

            Returns the couple: <specie_name>, <float>
        """

        regex_first_letter = r"[a-zA-Z]"
        matches = re.finditer(regex_first_letter, specie_in_reac)


        i_first_letter = -1
        for match in matches:
            i_first_letter = match.start()
            break

        if i_first_letter == -1:
            raise TypeError("No specie found in reactant/product {}.".format( specie_in_reac ) )

        specie_name = specie_in_reac[i_first_letter:]

        if specie_name not in self._species_name:
            raise TypeError("Reaction mentions specie {} which was not found in species list.".format( specie_name ) )

        specie_coeff_str = specie_in_reac[:i_first_letter]

        if specie_coeff_str == '':
            specie_coeff = 1
        else:
            specie_coeff = float( specie_coeff_str )

        return specie_name, specie_coeff

    def _get_reactions( self ):
        reactions = OrderedDict( {} )

        for reaction in self._json_dict[ 'reactions' ]:
            reactions[ reaction ] = self._parse_reaction( reaction )
        return reactions


class massActionCRN( object ):

    def __init__( self, name, species_name, x0, rates, reactions, debug = 0 ):
        """
            name: str
                Name of the CRN.
            species_name: list of str
                List of the name of species, the order matters as it fixes the
                order in all vectors of species (e.g. x0_np).
            x0: OrderedDict: str -> float
                Gives initial relative concentrations for each specie.
                The values should sum to one.
            rates: OrderedDict str -> float
                Gives rates for each reaction. Rates are named with str:
                    -k#i with i the index of the reaction
                    -k#-i if reverse rate.
            reactions: OrderedDict: str -> list of 2-tuples.
                Keys are reaction string e.g.: X+Y -> Z, and values are reaction list.
                These lists contain one or two element depending on `->` or `<->`.
                A reaction is a tuple of reactant/product relative stochiometry.
                e.g.: (np.array([-1,-1,0]),np.array([0,0,1]))
        """
        self._debug = debug
        self._debug_name = "[MA CRN]"

        self._species_name = species_name
        self._x0           = x0
        self._rates        = rates
        self._reactions    = reactions

        self._x0_np = np.array( list( self._x0.values() ) )
        self._check_x0_np()

        self._rates_np = np.array( list( self._rates.values() ) )

        self._reaction_list = self._get_reaction_list()
        self._reactants_matrix = np.vstack( [ r[ 0 ] for r in self._reaction_list] )
        self._gamma    = self._get_gamma()

        if self._debug > 0:
            print(self._debug_name+"species name:", self._species_name)
            print(self._debug_name+"x0:", self._x0)
            print(self._debug_name+"rates:", self._rates)
            print( self._debug_name+'reactions:')
            for r in self._reactions:
                print( r )
                for rr in self._reactions[ r ]:
                    print("\t", rr)

            print(self._debug_name+"x0_np:", self._x0_np)
            print(self._debug_name+"rates_np:", self._rates_np)
            print(self._debug_name+"reaction_list:", self._reaction_list)
            print(self._debug_name+"reactants_matrix:", self._reactants_matrix)
            print(self._debug_name+"gamma:", self._gamma)

    def _get_rate_list( self, i_reaction ):
        """ Return the rate list corresponding to reaction with index `i_reaction`.
        """
        rate_list = []
        rate_names = [ "k{}".format( i_reaction + 1 ),
                        "k-{}".format( i_reaction + 1 ) ]

        for rate_name in rate_names:
            if rate_name in self._rates:
                rate_list.append( self._rates[ rate_name ] )

        return rate_list

    def __str__( self ):
        """ Returns a pretty string summarizing the CRN.
        """
        to_return = "initial relative concentrations:\n"
        for specie in self._x0:
            to_return += "{}:{} ".format( specie, self._x0[specie] )

        to_return += "\n\nreactions:\n"

        for i, reaction in enumerate( self._reactions ):
            to_return += "{}. ".format(i+1) + reaction + " {}".format(self._get_rate_list( i ))

            if i + 1 != len( self._reactions ):
                to_return += "\n"

        #to_return += "\n\nrates:\n"
        #for i, rate_name in enumerate( self._rates ):
            #to_return += "{}: {}".format(rate_name,self._rates[rate_name])

            #if i + 1 != len( self._rates ):
                #to_return += "\n"

        return to_return

    def get_species_names( self ):
        return self._species_name

    def plot_dynamics( self, figsize = (20, 10) ):
        """ Plots the CRN dynamics.
        """
        history = self.integrate()
        plt.figure( figsize = figsize )

        for i, specie in enumerate( self.get_species_names() ):
            plt.plot( history[ :, i ], label = specie )

        plt.legend()
        plt.show()

    def build_UI( self, manual = True ):
        """ Builds an interactive UI for the CRN.
        """

        str_repr = str( self )
        system_label = widgets.Textarea( str_repr, disabled = True,
                                         rows = 1 + str_repr.count("\n") )

        def func( **kwargs ):
            global system_label

            new_crn = massActionCRN( self._species_name, self._x0, kwargs, self._reactions, self._debug )

            str_repr = str( new_crn  )
            system_label = widgets.Textarea( str_repr, disabled = True,
                                             rows = 1 + str_repr.count("\n") )

            new_crn.plot_dynamics()
            display( system_label )

        rates_widget = OrderedDict({})
        for rate in self._rates:
            widget = widgets.FloatSlider( value = self._rates[ rate ],
                                                min   = 1e-4,
                                                max   = 1e-2,
                                                step  = 1e-8,
                                                description = "rate {}".format( rate ),
                                                readout_format='.4f', )
            rates_widget[ rate ] = widget

        #w = interactive(func, {'manual': manual }, **rates_widget )
        output = interactive_output(func, rates_widget)
        output.layout.height = '720px'


        return widgets.VBox( [ output ] + list( rates_widget.values() ) )


    def integrate( self ):
        """ Integrates the dynamical system and gives the simulation results.
        """
        time = np.linspace(0, 20000, 10000)
        return integrate.odeint( self._mass_action_dynamic(), self._x0_np, time )

    def _mass_action_dynamic( self ):
        """ Implements the mass action dynamic. Returns a ready to use function
            for scipy to integrate.
        """
        def f(x, t0):
            speed_vector = (x**abs(self._reactants_matrix)).prod(axis = 1)*self._rates_np
            return self._gamma.dot( speed_vector )
        return f

    def _check_x0_np( self ):
        """ Check that self._x0_np is corresponding to relative concentrations.
        """
        if not np.isclose( self._x0_np.sum(), 1.0 ):
            raise TypeError("The CRN initial relative concentrations do not sum to 1: {}.".format(self._x0_np))

    def _get_reaction_list( self ):
        """ Transforms the reactions into a simpler object wich is a list
            of all the reactions.
        """
        return list( reduce( lambda a,b: a + b,
                             list( self._reactions.values() ) ) )

    def _get_gamma( self ):
        """
            Returns the stochiometric matrix associated to the CRN
        """
        return np.array( list( map( lambda x: reduce(lambda a,b: a+b, x),
                         self._reaction_list ) ) ).T
