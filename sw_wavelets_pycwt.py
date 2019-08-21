import glob
import os

def get_graph_from_file(in_filepath, out_folder, out_filename):
    # Get data
    p1, p2 = np.loadtxt(in_filepath, unpack=True)
    
    # TODO implement

    # Save graph to file
    # TODO implement
    plt.savefig('{}/{}.png'.format(out_folder, out_filename))
    # ----------------------------------------------
    # or show the graph
    #plt.show()

def get_graph_name_from_filepath(filepath):
    result = ""

    # Get only file name
    result = os.path.basename(filepath)

    # Truncate extension
    result = result[: result.rfind(".")]

    # Truncate all before station name
    result = result[result.find("wevelet_") + len("wevelet_") :]

    # Truncate passage number
    # Some tricky moment we need second from right _
    third_ = result.rfind("_")
    second_ = result.rfind("_", 0, third_) 
    result = result[:second_] + result[third_:]

    return result

def process_all_files_in_folder(in_folder, out_folder):
    for filepath in glob.glob(in_folder + "/" + "*.dat"):
        graph_name = get_graph_name_from_filepath(filepath)
        print(" >>> Process file '{}' with output name '{}'".format(filepath, graph_name))
        get_graph_from_file(filepath, out_folder, graph_name)

def main():
    process_all_files_in_folder("./input", "./output")

if __name__ == "__main__":
    main()