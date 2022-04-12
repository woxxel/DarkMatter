



# create class with simulator method for sbi
class sbi_darkmatter:

    def __init__():

        pass


    def register_data(self,pathToParameters):

        pass


    def create_prior(self):
        # uniform priors on all?
        # maybe rather Gaussian/student-t?

        """
            alpha: uniform
            tau_I (3x): student-t around empirical data
            tau_n: gauss around empirical?
            rate: student-t around measured stuff
        """

        pass


    def simulator(self,parameter_set):

        # should call c++ with as little overhead as possible (what can be avoided? what not?)

        pass
