from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
import subprocess

# Bounded region of parameter space
fcl_paramterbounds = {'StartFitSize': (3, 20), 
                      'NMissPoints':  (0, 15),
                      'MaxResidualDiff': (0,2),
                      'TrackMaxAdjacentSPDistance': (0,2).
                      'MaxAverageResidual': (0,4)
                      }


fcl_paramtertypes = {'StartFitSize': "int",
                     'NMissPoints': "int",
                     'MaxResidualDiff': "float",
                     'TrackMaxAdjacentSPDistance': "float",
                     'MaxAverageResidual': "float"
                     }

   
fcl_filename="SBNRecoValidation.fcl"
colon=":"
dot="."
newline="\n"

def black_box_function(StartFitSize, y):
    """Function with unknown internals we wish to maximize.

    This is just serving as an example, for all intents and
    purposes think of the internals of this function, i.e.: the process
    which generates its output values, as unknown.
    """
    return -StartFitSize ** 2 - (y - 1) ** 2 + 1


def replace_fcl_param(fcl_param, fcl_param_name):

    replace_fcl_param_name=dot+fcl_param_name+colon

    if(fcl_paramtertypes[fcl_param_name] == "int"):
        fcl_param=fcl_param.round()

    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(fcl_filename) as old_file:
            for line in old_file:
                new_fcl_param_name=line
                if(replace_fcl_param_name in line):
                    colon_loc = line.find(':')
                    new_fcl_param_name=line[:colon_loc+1]+str(fcl_param) + newline
                    line=line[colon_loc+1:]

                new_file.write(line.replace(line,new_fcl_param_name))

    #Remove original file
    remove(fcl_filename)
    #Move new file
    move(abs_path,fcl_filename)

optimizer = BayesianOptimization(
    f=None,
    pbounds=fcl_paramterbounds,
    verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=1,
)

utility = UtilityFunction(kind="ucb", kappa=2.5, xi=0.0)

for _ in range(5):

    #Get the suggested next parameters 
    next_point_to_probe = optimizer.suggest(utility)
    print("Next point to probe is:", next_point_to_probe)

    #Change the fcl config 
    for param in next_point_to_probe: 
        replace_fcl_param(next_point_to_probe.get(param),param)

    #Run the bash script required script must return the value to be maximised.
    bashCommand = "source ./larsoft_command.sh"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    target, error = process.communicate()
#    target = black_box_function(**next_point_to_probe)
    print("Found the target value to be:", target)
    
    #Add to the opttimsation processes 
    optimizer.register(
        params=next_point_to_probe,
        target=target,
    )

    print(target, next_point_to_probe)

print(optimizer.max)

