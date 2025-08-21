import cbor2
import numpy as np
import os
import gc
import itertools
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn.layout.venn3 import DefaultLayoutAlgorithm
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

class SearchNode:
    def __init__(self, bound_list, dur_list):
        self.bounds = set()
        self.times = list()
        self.children = []

        for b in bound_list:
            self.bounds.add(b)

        for dur in dur_list:
            nanos = dur['secs'] * 1e9 + dur['nanos']
            self.times.append(nanos)

    # Given cbor tree file, return root of SearchNode tree
    def from_cbor(tree_cbor):
        parent = SearchNode(tree_cbor['bounds'], tree_cbor['times'])

        for node in tree_cbor['children']:
            temp = SearchNode.from_cbor(node)
            parent.children.append(temp)

        return parent
    
    def bound_time(self, bound):
        index = 0
        match bound:
            case 'Log': index = 0
            case 'Int': index = 1
            case 'VecSimple': index = 2
            case 'VecSmallFrag': index = 3
            case 'Memoize': index = 4
        
        return self.times[index]
    
    # Given a list of bounds, tell the number of states searched
    # with these bounds enabled
    # Bounds can be Log, Int, VecSimple, VecSmallFrags
    def num_searched(self, bounds):
        if len(bounds & self.bounds) > 0:
            return 1
        else:
            return 1 + sum([x.num_searched(bounds) for x in self.children])
        
    # Prove score given bounds
    def score(self, bounds):
        total_weight = 500
        for b in bounds:
            total_weight += self.bound_time(b)
            if b in self.bounds:
                return total_weight
        
        return total_weight + sum([x.score(bounds) for x in self.children])

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
    pruned = dict()

    bounds = ['Int', 'VecSimple', 'VecSmallFrags']

    for mol in mols:
        print(mol)
        f = open(dataset + mol + '/tree.cbor', 'rb')
        data = cbor2.load(f)
        f.close()
        root = SearchNode.from_cbor(data)

        states_total = root.num_searched(set())
        total_search.append(states_total)

        for n in range(1, len(bounds) + 1):
            for comb in itertools.combinations(bounds, n):
                num_pruned = states_total - root.num_searched(set(comb))
                if comb in pruned:
                    pruned[comb].append(num_pruned)
                else:
                    pruned[comb] = [num_pruned]

        del data
        gc.collect()

    with open(outfile, 'wb') as f:
        cbor2.dump({'mols': mols, 'totals': total_search, 'pruned': pruned}, f)

'''
Broken for now
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
    plt.ylabel('# states pruned / # states searched (log-only)')
    plt.show()'''

def score_mol(mol_dir, bounds):
    exists = False

    # Try loading data if saved
    if os.path.exists(mol_dir + 'scores.cbor'):
        f = open(mol_dir + 'scores.cbor', 'rb')
        saved_bounds = cbor2.load(f)
        scores = cbor2.load(f)
        f.close()

        if set(bounds) <= set(saved_bounds):
            exists = True
    
    # Data was not saved
    if not exists:
        f = open(mol_dir + 'tree.cbor', 'rb')
        tree_data = cbor2.load(f)
        f.close()

        scores = []
        for n in range(1, len(bounds) + 1):
            for perm in itertools.permutations(bounds, n):
                score = score_cbor(tree_data, perm)
                scores.append([list(perm), score])
        
        scores.sort(key=lambda x: x[1])

        f = open(mol_dir + 'scores.cbor', 'wb')
        cbor2.dump(bounds, f)
        cbor2.dump(scores, f)
        f.close()

        # Explicitly call garbage collection to reduce mem usage
        del tree_data
        gc.collect()
    
    return scores

def score_cbor(node, bounds):
    total_weight = 0
    node_bounds = node['bounds']
    node_times = node['times']
    default_weight = node_times[5]['secs'] * 1e9 + node_times[5]['nanos']

    for b in bounds:
        index = 0
        match b:
            case 'Log': index = 0
            case 'Int': index = 1
            case 'VecSimple': index = 2
            case 'VecSmallFrags': index = 3
            case 'Memoize': index = 4

        bound_dur = node_times[index]
        bound_time = bound_dur['secs'] * 1e9 + bound_dur['nanos']
        
        total_weight += bound_time
        if b in node_bounds:
            return total_weight

    return total_weight + default_weight + sum([score_cbor(x, bounds) for x in node['children']])

def venn_diagram(inputfile):
    # input file is generated by output_states_searched
    f = open(inputfile, 'rb')
    data = cbor2.load(f)
    f.close()

    mols = data['mols']
    total_search = data['totals']
    pruned = data['pruned']

    num_mols = len(mols)
    # Create matrix to go from states_searched output to 
    # venn diagram values
    include_mat = np.array([
        [1, 0, 1, 0, 1, 0, 1],
        [0, 1, 1, 0, 0, 1, 1],
        [1, 1, 1, 0, 1, 1, 1],
        [0, 0, 0, 1, 1, 1, 1],
        [1, 0, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1],
    ])
    include_mat = np.linalg.inv(include_mat)

    circles = np.zeros(7)

    for i in range(num_mols):
        #Abc, aBc, ABc, abC, AbC, aBC, ABC
        pruned_list = np.zeros(7)
        pruned_list[0] = pruned[('Int',)][i]
        pruned_list[1] = pruned[('VecSimple',)][i]
        pruned_list[2] = pruned[('Int', 'VecSimple',)][i]
        pruned_list[3] = pruned[('VecSmallFrags',)][i]
        pruned_list[4] = pruned[('Int', 'VecSmallFrags',)][i]
        pruned_list[5] = pruned[('VecSimple', 'VecSmallFrags',)][i]
        pruned_list[6] = pruned[('Int', 'VecSimple', 'VecSmallFrags',)][i]

        # Change representation
        pruned_list = include_mat @ pruned_list
        
        circles = np.add(circles, pruned_list)

    venn3(subsets = circles, 
          set_labels=('Int', 'VecSimple', 'VecSmallFrags'), 
          layout_algorithm=DefaultLayoutAlgorithm(fixed_subset_sizes=(1,1,1,1,1,1,1)))
    plt.show()

def gen_tree(tree, child_names, parent_name):
    for cname in child_names:
        name = parent_name + [cname]

        # Create child node
        child = tree.add_child(name=name)

        # Add name to node
        child.add_face(TextFace(cname), column=0, position='branch-top')

        # Recurse
        new_list = child_names.copy()
        new_list.remove(cname)
        gen_tree(child, new_list, name)

def bound_order_tree(dir, bounds, dataset=False):
    # Generate tree structure
    t = Tree(name='root')
    t.add_face(TextFace('None'), column=0, position='branch-top')
    ts = TreeStyle()
    ts.show_scale = False
    ts.scale = 240
    ts.branch_vertical_margin = 30
    ts.show_leaf_name = False
    gen_tree(t, bounds, list())

    # Load scores
    scores = list()
    if dataset == True:
        mols = os.listdir(dir)

        totals = []
        for n in range(1, len(bounds) + 1):
            for perm in itertools.permutations(bounds, n):
                totals.append([list(perm), 0])

        for mol in mols:
            print(mol)
            mol_scores = score_mol(dir + mol + '/', bounds)

            low = mol_scores[0][1]
            high = mol_scores[-1][1]

            for score in mol_scores:
                score_name = score[0]
                score_val = score[1]

                prop_val = (score_val - low) / (high - low)

                i = 0
                while totals[i][0] != score_name:
                    i += 1
                
                totals[i][1] += score_val

        totals.sort(key=lambda x: x[1])
        scores = totals
    else:
        scores = score_mol(dir, bounds)

    low = scores[0][1]
    high = scores[-2][1]

    # Add visual elements to tree
    for n in t.traverse():
        nstyle = NodeStyle()
        color = '#000000'
        default_rbg = (0, 128, 255)
        name = n.name

        if name != 'root':
            # Find score entry for this name
            i = 0
            while scores[i][0] != n.name:
                i += 1
            
            # Find proportional score
            val = scores[i][1]
            prop_val = 1 - (val - low) / (high - low)

            # Add score to node
            n.add_face(TextFace(round(prop_val, 3)), column=0, position='branch-bottom')

            # Color according to proporional score
            color = '#'
            for i in default_rbg:
                round_val = round(prop_val * i)
                color += f'{round_val:02X}'

        nstyle['shape'] = 'square'
        nstyle['size'] = 60
        nstyle['fgcolor'] = color

        n.set_style(nstyle)

    t.show(tree_style=ts)


bound_order_tree('eval_out/coconut_55/', ['Log', 'Int', 'VecSimple', 'VecSmallFrags'], dataset=True)

'''
low = (0, 0, 255)
mlow = (0, 255, 0)
mhigh = (255, 255, 0),
high = (255, 0, 0)
'''