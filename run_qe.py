from os import system

def run_pwx(pwin,pwout):
    """
    Runs pw.x serially on my laptop
    To be modified for a different machine
    """
    cmd = "/Users/dariorocca/work/qe-6.3/bin/pw.x < "+pwin+" > "+pwout
    system(cmd)
