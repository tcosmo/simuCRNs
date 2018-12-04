import json
import simuCRNs.massActionCRN
import simuCRNs.uiCRNs

from collections import OrderedDict

def fromJSON( json_file_name, debug = 0 ):
    """ Constructs a CRN object from a JSON specification.

        debug: depth of required debugging output. If 0, no debug output.

        Example for the mass action model:

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
    """

    with open( json_file_name, "r" ) as f:
        json_dict = json.load( f,  object_pairs_hook=OrderedDict )

    if not "type" in json_dict:
        raise TypeError("No simulation type specified: should be \
                        \"type\": \"mass action\" or \"type\": \"stochastic\".")

    if json_dict["type"] == "mass action":
        return simuCRNs.massActionCRN.massActionJSON( json_dict, debug ).parse()
    elif json_dict["type"] == "stochastic":
        raise NotImplementedError("The stochastic model is not yet implemented.")
    else:
        raise TypeError("CRN type {} not known: should be \
                        \"type\": \"mass action\" or \"type\": \"stochastic\".".format(json_dict["type"]))
