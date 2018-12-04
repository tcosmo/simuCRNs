import simuCRNs.massActionCRN

from ipywidgets import interact, interactive, fixed, interact_manual, interactive_output
import scipy.integrate as integrate
import ipywidgets as widgets

from collections import OrderedDict

class defaultMassActionUI( object ):
    """ Default UI for Mass Action CRNs.
    """
    #TODO: add time management

    def __init__( self, mass_action_crn, default_sliders_config, sliders_config = {},
                        step_for_float_text = 0.1, ui_horizontal = True,
                        figsize = (8, 6), dpi = 80, layout_height = '740px' ):
            """ Specifies UI parameters.

                mass_action_crn: massActionCRN

                default_sliders_config: tuple
                    Default (min,max,step) for sliders.
                sliders_config: dict: str -> tuple
                    Specifies the (min, max, step) of slider of rate `rate_name`.
                step_for_float_text: float
                    Increment step for x0 `FloatText` widgets.
                ui_horizontal: bool
                    If True then the controlers will be displayed horizontally to
                    the plot. If False, vertically.
                figsize: tuple
                    Figure size for `plot_dynamics`.
                dpi: int
                    Figure dpi for `plot_dynamics`.
                layout_height: str
                    Height of the plot layout, useful to adjust for small / big screens.
            """

            self._crn = mass_action_crn
            self._default_sliders_config = default_sliders_config
            self._sliders_config = sliders_config
            self._step_for_float_text = step_for_float_text
            self._ui_horizontal = ui_horizontal
            self._figsize = figsize
            self._dpi = dpi
            self._layout_height = layout_height


    def build_UI( self ):
        """ Builds default UI for mass action CRNs.
        """
        str_repr = str( self )
        system_label = widgets.Textarea( str_repr, disabled = True,
                                         rows = 1 + str_repr.count("\n") )

        def func( **kwargs ):
            global system_label

            new_rates = OrderedDict( {} )
            for rate in self._crn.get_rates():
                new_rates[ rate ] = kwargs[ rate ]

            new_x0 = OrderedDict( {} )
            for specie in self._crn.get_species_names():
                new_x0[ specie ] = kwargs[ '#'+specie ]

            new_crn = simuCRNs.massActionCRN.massActionCRN( self._crn.get_name(), self._crn.get_species_names(), new_x0,
                                     new_rates, self._crn.get_reactions(), self._crn._debug )

            str_repr = str( new_crn  )
            system_label = widgets.Textarea( str_repr, disabled = True,
                                             rows = 1 + str_repr.count("\n") )

            new_crn.plot_dynamics( self._figsize, self._dpi )
            display( system_label )

        rates_widget = OrderedDict({})
        rates_widget_unit = OrderedDict({})
        for rate in self._crn.get_rates():
            config = self._default_sliders_config
            if rate in self._sliders_config:
                config = self._sliders_config[ rate ]

            widget = widgets.FloatSlider( value = self._crn.get_rates()[ rate ],
                                          min   = config[ 0 ],
                                          max   = config[ 1 ],
                                          step  = config[ 2 ],
                                          description = "rate {}".format( rate ),
                                          readout_format='.2e', )
            rates_widget[ rate ] = widget
            rates_widget_unit[ rate ] = widgets.HBox( [ widget,
                widgets.HTMLMath(
                    value=self._crn.get_rate_unit( rate ),
                    ) ] )

        x0_widget = OrderedDict( {} )
        for specie in self._crn.get_species_names():
            widget = widgets.BoundedFloatText( value = self._crn.get_x0()[ specie ],
                                               min = 0.0,
                                               step = self._step_for_float_text,
                                               description = specie + " (Mol)" )
            x0_widget[ '#'+specie ] = widget




        all_widgets = OrderedDict(list(rates_widget.items()) +
                                  list(x0_widget.items()))
        all_widgets_units = OrderedDict(list(rates_widget_unit.items()) +
                                  list(x0_widget.items()))
        #w = interactive(func, {'manual': manual }, **rates_widget )
        output = interactive_output(func, all_widgets)
        output.layout.height = self._layout_height


        encapsulation = widgets.HBox if self._ui_horizontal else widgets.VBox

        return encapsulation( [ output ] +
                              [ widgets.VBox( list( all_widgets_units.values() ) ) ] )
