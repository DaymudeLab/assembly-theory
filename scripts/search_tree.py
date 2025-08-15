import cbor2
import os
import matplotlib.pyplot as plt
import itertools

class SearchNode:
    def __init__(self, bound_list):
        self.bounds = set()
        self.children = []

        for b in bound_list:
            self.bounds.add(b)

    # Given cbor tree file, return root of SearchNode tree
    def from_cbor(tree_cbor):
        parent = SearchNode(tree_cbor['bounds'])

        for node in tree_cbor['children']:
            temp = SearchNode.from_cbor(node)
            parent.children.append(temp)

        return parent
    
    # Given a list of bounds, tell the number of states searched
    # with these bounds enabled
    # Bounds can be Log, Int, VecSimple, VecSmallFrags
    def num_searched(self, bounds):
        if len(bounds & self.bounds) > 0:
            return 1
        else:
            return 1 + sum([x.num_searched(bounds) for x in self.children])
        
    # Prove score given bounds and their timing information
    # bounds_weights: [(bound, weight)]
    def score(self, bounds_weights):
        total_weight = 0
        for b, w in bounds_weights:
            total_weight += w
            if b in self.bounds:
                return total_weight
        
        return total_weight + sum([x.score(bounds_weights) for x in self.children])

class Timer:
    def __init__(self, timer, bound_name):
        # timer: {metric -> [{time_unit -> time}, num]}
        self.timer = timer
        self.bound = bound_name

    def from_cbor(name, timer_data):
        name_map = {
            'Log': 'log_timer',
            'Int': 'int_timer',
            'VecSimple': 'vec_simple_timer',
            'VecSmallFrags': 'vec_small_frags_timer',
            'Memoize': 'memoize_timer',
        }

        bound_name = name_map[name]
        return Timer(timer_data[bound_name], name)

    # Average time of the timer in nanos
    def average(self):
        total_time = 0
        total_num = 0
        for val in self.timer.values():
            total_time += 1e9 * val[0]['secs'] + val[0]['nanos']
            total_num += val[1]

        return total_time / total_num

def print_average_bound_times(dataset):
    timer_strings = ['Log', 'Int', 'VecSimple', 'VecSmallFrags']

    total_avg = {name: 0 for name in timer_strings}
    num_mols = 0

    for mol_name in os.listdir(dataset):
        num_mols += 1

        timer_file = open(dataset + mol_name + '/timer.cbor', 'rb')
        timers = cbor2.load(timer_file)
        timer_file.close()

        for name in timer_strings:
            timer = Timer.from_cbor(name, timers)
            avg = timer.average()
            total_avg[name] += avg

    for key, val in total_avg.items():
        print(key, ":", val / num_mols)

def output_states_searched(dataset, outfile):
    # list of mols to gather data for
    mols = os.listdir(dataset)
    # Total states searched for each mol
    total_search = []
    # Total pruned for each bound
    # Each bound also includes all previous
    ordered_bounds = {
        'Int': [],
        'VecSimlpe': [],
        'VecSmallFrags': [],
    }

    for mol in mols:
        print(mol)
        f = open(dataset + mol + '/tree.cbor', 'rb')
        data = cbor2.load(f)
        f.close()
        root = SearchNode.from_cbor(data)

        bounds_set = set()
        states_total = root.num_searched(bounds_set)
        total_search.append(states_total)
        for bound in ordered_bounds:
            bounds_set.add(bound)
            ordered_bounds[bound].append(states_total - root.num_searched(bounds_set))

    with open(outfile, 'wb') as f:
        cbor2.dump({'mols': mols, 'totals': total_search, 'pruned': ordered_bounds}, f)

def plot_states_searched(infile):
    f = open(infile, 'rb')
    data = cbor2.load(f)
    f.close()

    mols = data['mols']
    total_search = data['totals']
    ordered_bounds = data['pruned']

    # The list of lists to be plotted as a stacked bar chart
    # Modify this to change what is plotted
    to_plot = dict()
    for bound, pruned in ordered_bounds.items():
        to_plot[bound] = [x / y for x, y in zip(pruned, total_search)]

    big_list = []
    for i in range(len(mols)):
        small_list = []
        for _, pruned in to_plot.items():
            small_list.append(pruned[i])
        big_list.append(small_list)

    big_list.sort(key=lambda x: x[len(to_plot) - 1])
    for (i,bound) in enumerate(to_plot.keys()):
        to_plot[bound] = [l[i] for l in big_list]


    _, ax = plt.subplots()
    bottom = [0] * len(mols)
    for label, y_vars in to_plot.items():
        top = [x[0] - x[1] for x in zip(y_vars, bottom)]
        ax.bar(mols, top, label=label, bottom = bottom)
        bottom = y_vars

    plt.xticks([])
    plt.legend(loc="upper left")
    plt.xlabel('Molecule')
    plt.ylabel('% States Pruned')
    plt.show()

def score_mol(mol_file, bounds):
    f = open(mol_file + '/timer.cbor', 'rb')
    timer_data = cbor2.load(f)
    f.close()
    timers = [Timer.from_cbor(b, timer_data) for b in bounds]

    f = open(mol_file + '/tree.cbor', 'rb')
    tree_data = cbor2.load(f)
    f.close()
    root = SearchNode.from_cbor(tree_data)

    bounds_weights = [(x.bound, x.average()) for x in timers]
    scores = []
    for n in range(1, len(bounds_weights) + 1):
        for perm in itertools.permutations(bounds_weights, n):
            bounds = [x[0] for x in perm]
            score = root.score(perm)
            scores.append((bounds, score))
    
    scores.sort(key=lambda x: -x[1])
    for s in scores:
        print(s)


dataset = 'eval_out/coconut_55/'
mol = dataset + '16.mol'
score_mol(mol, ['Log', 'Int', 'VecSimple', 'VecSmallFrags'])


