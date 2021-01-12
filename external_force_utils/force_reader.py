import forces

def read_force_file(file):
    force_list = []
    with open(file, 'r') as f:
        n = 0
        l = f.readline()
        while l:
            if l == "{{":  #new force
                l = f.readline()
                args = {}
                while l != "}}": #read until end of description
                    l = l.split("=")
                    if l[0].strip() == "type":
                        t = l[0]
                    else:
                        args[l[0].strip()] = l[1].strip()
                    l = f.readline()
                force_list.append(getattr(forces, t)(**args)) #calls the function "t" from module "forces"
            l = f.readline()
    return(force_list)

"""
Write a list of forces out to a file.

Parameters:
    force_list (list): a list of force dictionaries
    filename (str): file to write out to
    <optional> mode (str): the mode of the python open funciton.  Defaults to 'w'
    change mode to 'w+' if you want to append instead of overwrite.
"""
def write_force_file(force_list, filename, mode='w'):
    with open(filename, mode=mode) as f:
        for force in force_list:
            out = "{{\n"
            for k in force.keys():
                out += "{} = ".format(k)
                if isinstance(force[k], list):
                    out += ", ".join(force[k])
                else:
                    out += force[k]
                out += "\n"
            out += "}}\n\n"
        f.write(out)