from pySDC.tutorial.step_6.A_run_non_MPI_controller import main as main_A


def main():
    """
    A simple test program to do check PFASST for odd numbers of processes
    """
    main_A(num_proc_list=[3, 5, 7, 9], fname='step_6_B_out.txt', multi_level=True)


if __name__ == "__main__":
    main()
