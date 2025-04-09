"""
The network size and the travel demand have been scaled, scale size may be:
1. the capacity and the travel demand is divided by 1000.
2. the free flow time is divided by 100.
The result still fails to align with the result in the paper.
However, through checking for many times, we believe that the logic is right.
Factors leading to the difference may be the scale size.
"""


import re
import time
import numpy as np
import matplotlib.pyplot as plt
from math import inf
from matplotlib.ticker import PercentFormatter


class Node:
    def __init__(self, node_id):
        self.node_id: int = node_id
        self.upstream_link: list[Link] = list()
        self.downstream_link: list[Link] = list()
        self.parent: Node = self
        self.dist: float = inf

    def __repr__(self):
        return f'NODE {self.node_id}'


class Link:
    def __init__(self, link_id, tail=None, head=None, capacity=None,
                 length=None, free_flow_time=None, alpha=None, beta=None):
        self.link_id: int = link_id
        self.tail: Node = tail
        self.head: Node = head
        self.capacity: float = capacity
        self.length: float = length
        self.fft: float = free_flow_time
        self.alpha: float = alpha
        self.beta: float = beta
        self.ue_flow: float = 0
        self.auxiliary_ue_flow: float = 0
        self.so_flow: float = 0
        self.auxiliary_so_flow: float = 0
        self.cost: float = 0
        self.marginal_cost: float = 0

    def __repr__(self):
        return (f'LINK {self.link_id} cost = {self.cost}, marginal_cost = {self.marginal_cost},'
                f' so_flow = {self.so_flow}, ue_flow = {self.ue_flow}')

    def update_cost(self) -> None:
        self.cost = self.fft * (1 + self.alpha * ((self.ue_flow + self.so_flow) / self.capacity) ** self.beta)

    def update_marginal_cost(self) -> None:
        self.marginal_cost = (self.fft * (1 + self.alpha * ((self.ue_flow + self.so_flow) / self.capacity) ** self.beta)
                              + self.fft * self.alpha * self.beta * (
                                      (self.ue_flow + self.so_flow) / self.capacity) ** self.beta)

    def get_cost(self, ue_param, so_param) -> float:
        return self.fft * (1 + self.alpha * ((ue_param + so_param) / self.capacity) ** self.beta)

    def get_marginal_cost(self, ue_param, so_param) -> float:
        return (self.fft * (1 + self.alpha * ((ue_param + so_param) / self.capacity) ** self.beta)
                + self.fft * self.alpha * self.beta * ((ue_param + so_param) / self.capacity) ** self.beta)


class ODPair:
    def __init__(self, origin, destination, demand):
        self.origin: Node = origin
        self.destination: Node = destination
        self.total_demand: float = demand
        self.ue_demand: float = 0
        self.so_demand: float = 0

    def __repr__(self):
        return f'ODPair {self.origin.node_id}->{self.destination.node_id}={self.total_demand}'


class Network:
    def __init__(self, name, sst=1, od_scale=1e-3, net_scale=1e-2):
        self.name: str = name
        self.sst = sst  # demand sensitivity
        self.od_scale = od_scale
        self.net_scale = net_scale
        self.Nodes: list[Node] = list()
        self.Links: list[Link] = list()
        self.ODPairs: list[ODPair] = list()
        self.num_node: int = 0
        self.num_link: int = 0
        self.total_flow: float = 0
        self.read_network()
        self.read_OD()

    def read_network(self):
        path = f'{self.name}\\{self.name}_net.txt'
        with open(path, 'r', encoding='UTF-8') as file:
            lines = file.readlines()
        pattern = re.compile(r'[\w.~]+')
        data = [pattern.findall(line) for line in lines if len(pattern.findall(line)) != 0]
        self.num_node, self.num_link = int(data[1][-1]), int(data[3][-1])
        for i in range(len(data)):
            if '~' in data[i] and "ORIGINAL" not in data[i]:
                data = data[i + 1:]
                break
        self.Nodes = [Node(i) for i in range(self.num_node + 1)]
        self.Links = [Link(0)]
        for index, line in enumerate(data):
            cur_link = Link(index+1, self.Nodes[int(line[0])], self.Nodes[int(line[1])], float(line[2])*self.od_scale,
                            float(line[3]), float(line[4])*self.net_scale, float(line[5]), float(line[6]))
            self.Links.append(cur_link)
            self.Nodes[int(line[0])].downstream_link.append(cur_link)
            self.Nodes[int(line[1])].upstream_link.append(cur_link)

    def read_OD(self):
        path = f'{self.name}\\{self.name}_trips.txt'
        with open(path, 'r', encoding='UTF-8') as file:
            # Process the text file
            lines = file.readlines()
        pattern = re.compile(r'[0-9.]+|Origin')
        data = [pattern.findall(line) for line in lines if len(pattern.findall(line)) != 0]
        self.total_flow = float(data[1][0]) * self.od_scale
        for i in range(len(data)):
            if 'Origin' in data[i]:
                data = data[i:]
                break
        origin = None
        for line in data:
            if "Origin" in line:
                origin = self.Nodes[int(line[-1])]
                continue
            for i in range(len(line) // 2):
                destination = self.Nodes[int(line[2 * i])]
                demand = float(line[2 * i + 1]) * self.sst * self.od_scale
                if demand != 0:
                    self.ODPairs.append(ODPair(origin, destination, demand))
        if abs(self.total_flow - sum([od.total_demand for od in self.ODPairs])) > 1:
            raise ValueError("Demand in the file doesn't match with the total OD flow.")


# Solve the lower problem
class MixedEquilibrium:
    def __init__(self, network, ue_demand, so_demand, ue_gap, so_gap, FW_max_iter, LS_max_iter):
        self.network: Network = network
        self.ue_demand: list[float] = ue_demand
        self.so_demand: list[float] = so_demand
        self.ue_gap: float = ue_gap
        self.so_gap: float = so_gap
        self.FW_max_iter: int = FW_max_iter  # FW algorithm convergence gap
        self.LS_max_iter: int = LS_max_iter  # Line search convergence gap
        self.TSTT: float = 0

    # Main method, return the total system travel time
    def run(self) -> tuple[float, np.ndarray, np.ndarray]:
        # Initialize
        iteration, cur_ue_gap, cur_so_gap = 0, inf, inf
        self.initialize()
        while iteration < self.FW_max_iter and (cur_ue_gap > self.ue_gap or cur_so_gap > self.so_gap):
            # All or nothing assignment for UE users
            ue_search_dir = self.all_or_nothing_ue()
            # All or nothing assignment for SO users
            so_search_dir = self.all_or_nothing_so()
            # Check convergence for every 200 iterations to improve efficiency
            if iteration % 200 == 0:
                cur_ue_gap = self.check_ue_convergence()
                cur_so_gap = self.check_so_convergence()
            # Line search
            optimal_step = self.line_search(ue_search_dir, so_search_dir)
            for i, link in enumerate(self.network.Links[1:]):
                link.ue_flow += optimal_step * ue_search_dir[i]
                link.so_flow += optimal_step * so_search_dir[i]
            iteration += 1

        total_travel_time = sum([(link.cost * (link.so_flow + link.ue_flow)) for link in self.network.Links[1:]])
        ue_flow_pattern = np.array([link.ue_flow for link in self.network.Links[1:]])
        so_flow_pattern = np.array([link.so_flow for link in self.network.Links[1:]])
        return total_travel_time, ue_flow_pattern, so_flow_pattern

    # Load the SO and UE flow; Initialize cost and marginal cost; Initialize the link flow based on all-or-nothing.
    def initialize(self):
        for index, od in enumerate(self.network.ODPairs):
            od.ue_demand = self.ue_demand[index]
            od.so_demand = self.so_demand[index]
        for link in self.network.Links[1:]:
            link.ue_flow, link.so_flow = 0, 0
            link.update_cost()
            link.update_marginal_cost()
        self.all_or_nothing_so()
        self.all_or_nothing_ue()
        for link in self.network.Links[1:]:
            link.ue_flow = link.auxiliary_ue_flow
            link.so_flow = link.auxiliary_so_flow

    def all_or_nothing_ue(self) -> np.ndarray:
        for link in self.network.Links[1:]:
            link.update_cost()
            link.auxiliary_ue_flow = 0
        for od in self.network.ODPairs:
            origin, destination, ue_demand = od.origin.node_id, od.destination.node_id, od.ue_demand
            shortest_path = self.dijkstra(origin, destination, marginal=False)
            for link in shortest_path:
                link.auxiliary_ue_flow += ue_demand
        return np.array([(link.auxiliary_ue_flow - link.ue_flow) for link in self.network.Links[1:]])

    def all_or_nothing_so(self) -> np.ndarray:
        for link in self.network.Links[1:]:
            link.update_marginal_cost()
            link.auxiliary_so_flow = 0
        for od in self.network.ODPairs:
            origin, destination, so_demand = od.origin.node_id, od.destination.node_id, od.so_demand
            shortest_path = self.dijkstra(origin, destination, marginal=True)
            for link in shortest_path:
                link.auxiliary_so_flow += so_demand
        return np.array([(link.auxiliary_so_flow - link.so_flow) for link in self.network.Links[1:]])

    def check_ue_convergence(self) -> float:
        if sum([od.ue_demand for od in self.network.ODPairs]) == 0:
            return -1
        SPTT = 0
        for od in self.network.ODPairs:
            min_path = self.dijkstra(od.origin.node_id, od.destination.node_id, marginal=False)
            min_dist = sum([link.cost for link in min_path])
            SPTT += min_dist * od.ue_demand
        TSTT = sum([link.ue_flow * link.cost for link in self.network.Links[1:]])
        return (TSTT / SPTT) - 1

    def check_so_convergence(self) -> float:
        if sum([od.so_demand for od in self.network.ODPairs]) == 0:
            return -1
        SPTT = 0
        for od in self.network.ODPairs:
            min_path = self.dijkstra(od.origin.node_id, od.destination.node_id, marginal=True)
            min_dist = sum([link.marginal_cost for link in min_path])
            SPTT += min_dist * od.so_demand
        TSTT = sum([link.so_flow * link.marginal_cost for link in self.network.Links[1:]])
        return (TSTT / SPTT) - 1

    def line_search(self, ue_dir, so_dir) -> float:
        iteration_ls, left, right, mid = 0, 0, 1, 0.5
        ue_flow_vector = np.array([link.ue_flow for link in self.network.Links[1:]])
        so_flow_vector = np.array([link.so_flow for link in self.network.Links[1:]])
        ue_flow_alpha = ue_flow_vector + mid * ue_dir
        so_flow_alpha = so_flow_vector + mid * so_dir
        cost_alpha = np.array([link.get_cost(ue_flow_alpha[i], so_flow_alpha[i])
                               for i, link in enumerate(self.network.Links[1:])])
        marginal_cost_alpha = np.array([link.get_marginal_cost(ue_flow_alpha[i], so_flow_alpha[i])
                                        for i, link in enumerate(self.network.Links[1:])])
        sigma = np.dot(cost_alpha, ue_dir) + np.dot(marginal_cost_alpha, so_dir)
        while (iteration_ls < self.LS_max_iter and (right - left) > self.so_gap) or sigma > 0:
            if sigma < 0:
                left = mid
            else:
                right = mid
            mid = 0.5 * (right + left)
            ue_flow_alpha = ue_flow_vector + mid * ue_dir
            so_flow_alpha = so_flow_vector + mid * so_dir
            cost_alpha = np.array([link.get_cost(ue_flow_alpha[i], so_flow_alpha[i])
                                   for i, link in enumerate(self.network.Links[1:])])
            marginal_cost_alpha = np.array([link.get_marginal_cost(ue_flow_alpha[i], so_flow_alpha[i])
                                            for i, link in enumerate(self.network.Links[1:])])
            sigma = np.dot(cost_alpha, ue_dir) + np.dot(marginal_cost_alpha, so_dir)
            iteration_ls += 1
        return mid

    """
    Find the shortest path based on the cost or the marginal cost.
    Set `marginal=False` for the former and `marginal=True` for the latter.
    """

    def dijkstra(self, o_id: int, d_id: int, marginal: bool) -> list[Link]:
        # initialize
        for node in self.network.Nodes:
            node.parent = node
            node.dist = inf
        self.network.Nodes[o_id].parent = -1
        self.network.Nodes[o_id].dist = 0
        # main loop
        SEL = [self.network.Nodes[o_id]]
        while SEL:
            SEL.sort(key=lambda n: n.dist, reverse=True)
            cur = SEL.pop()
            if cur.node_id == d_id:
                break
            for link in cur.downstream_link:
                if marginal:
                    potential_dist = cur.dist + link.marginal_cost
                else:
                    potential_dist = cur.dist + link.cost
                if link.head.dist > potential_dist:
                    link.head.parent = cur
                    link.head.dist = potential_dist
                    if link.head not in SEL:
                        SEL.append(link.head)
        # obtain the shortest path
        shortest_path = []
        cur_node = self.network.Nodes[d_id]
        while cur_node.parent != -1:
            p = cur_node.parent
            for link in cur_node.upstream_link:
                if link.tail == p:
                    temp = link
                    break
            else:
                raise ValueError("No shortest path is found.")
            shortest_path.append(temp)
            cur_node = p
        shortest_path.reverse()
        return shortest_path


# Solve the upper problem
class ORCS:
    def __init__(self, network, control_intensity, penalty_param, control_potential, gap1, gap2, gap3, gap4, max_iter1,
                 max_iter2, max_iter3, max_iter4):
        self.network: Network = network
        self.control_intensity: float = control_intensity
        self.penalty_param: float = penalty_param
        self.control_potential: float = control_potential
        self.gap1: float = gap1  # UE gap for lower problem
        self.gap2: float = gap2  # SO gap for upper problem
        self.gap3: float = gap3  # main gap for ORCS
        self.gap4: float = gap4  # for ADMM
        self.max_iter1: int = max_iter1  # Lower problem main loop max iteration
        self.max_iter2: int = max_iter2  # Lower problem line search max iteration
        self.max_iter3: int = max_iter3  # Upper problem main loop max iteration
        self.max_iter4: int = max_iter4  # Upper problem ADMM max iteration
        self.tstt_list: list = list()
        self.so_percent_list: list = list()

    def run(self):
        # initialize
        iteration, cur_gap3, num_od = 0, inf, len(self.network.ODPairs)
        shifted_flow = np.zeros(num_od)
        cur_ue_demand = np.array([od.total_demand for od in self.network.ODPairs])
        cur_so_demand = np.zeros(num_od)
        lp = MixedEquilibrium(self.network, cur_ue_demand, cur_so_demand, self.gap1, self.gap2, self.max_iter1,
                              self.max_iter2)
        lp.run()
        # main loop
        while iteration < self.max_iter3 and cur_gap3 > self.gap3:
            # calculate the approximation of gradient
            for link in self.network.Links[1:]:
                link.update_cost()
                link.update_marginal_cost()
            gradient = list()
            for od in self.network.ODPairs:
                origin, destination = od.origin.node_id, od.destination.node_id
                ue_mc = sum([link.marginal_cost for link in lp.dijkstra(origin, destination, marginal=False)])
                so_mc = sum([link.marginal_cost for link in lp.dijkstra(origin, destination, marginal=True)])
                gradient.append(so_mc-ue_mc)
            gradient = np.array(gradient)
            # solve the optimal control problem using ADMM
            ADMM_iteration, cur_gap1, cur_gap2 = 0, inf, inf
            c, u = np.zeros(num_od), np.zeros(num_od)
            new_shifted_flow = c - u - gradient / self.penalty_param
            while ADMM_iteration < self.max_iter4 and (cur_gap1 > self.gap4 or cur_gap2 > self.gap4):
                # update q
                new_shifted_flow = c - u - gradient / self.penalty_param
                # update c
                new_c = list()
                for i, od in enumerate(self.network.ODPairs):
                    temp1 = new_shifted_flow[i] + u[i]
                    temp2 = self.control_intensity / self.penalty_param
                    max_flow_shifted, min_flow_shifted = cur_ue_demand[i], -cur_so_demand[i]
                    if temp1 > temp2:
                        new_c.append(min(float(max_flow_shifted), float(temp1-temp2)))
                    elif -temp2 <= temp1 <= temp2:
                        new_c.append(0)
                    else:
                        new_c.append(max(float(min_flow_shifted), float(temp1+temp2)))
                new_c = np.array(new_c)
                # update u
                new_u = u + new_shifted_flow - new_c
                cur_gap1 = np.linalg.norm(new_shifted_flow-new_c, ord=1)
                cur_gap2 = np.linalg.norm(new_c-c, ord=1)
                c, u = new_c, new_u
                ADMM_iteration += 1
            cur_ue_demand -= new_shifted_flow
            cur_so_demand += new_shifted_flow
            # solve the new mixed traffic equilibrium
            lp = MixedEquilibrium(self.network, cur_ue_demand, cur_so_demand, self.gap1, self.gap2, self.max_iter1,
                                  self.max_iter2)
            tstt, ue_flow, so_flow = lp.run()
            cur_gap3 = np.linalg.norm(new_shifted_flow, ord=1)
            shifted_flow += new_shifted_flow
            iteration += 1
            so_percent = sum(shifted_flow) / self.network.total_flow
            self.tstt_list.append(tstt)
            self.so_percent_list.append(so_percent)
            print(f'Iteration {iteration}: cur_gap3 = {cur_gap3:.5f},'
                  f' TSTT = {tstt:.2f}, controlled ratio = {so_percent:.3f},'
                  f' total controlled demand = {sum(shifted_flow):.3f}')

    def plot_basic_data(self):
        ue_tstt, so_tstt = 74.80225, 71.94922
        num_iter = len(self.tstt_list)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        axes[0].plot(range(1, num_iter + 1), self.tstt_list)
        axes[0].set_ylabel("Total travel time", fontsize=12)
        axes[0].set_xlabel("Number of iteration", fontsize=12)
        axes[0].set_xticks(range(1, num_iter + 1))
        axes[0].set_ylim(71.8, 75.0)
        axes[0].set_xticklabels([str(i) for i in range(1, num_iter + 1)])
        axes[0].hlines(ue_tstt, xmax=num_iter, xmin=1, color="#FFA500", linestyles="--", label="UE")
        axes[0].hlines(so_tstt, xmax=num_iter, xmin=1, color="green", linestyles="--", label="SO")
        axes[0].legend(bbox_to_anchor=(0.97, 0.9))
        axes[1].plot(range(1, num_iter + 1), self.so_percent_list)
        axes[1].set_ylabel("% of SO users", fontsize=12)
        axes[1].set_xlabel("Number of iteration", fontsize=12)
        axes[1].set_xticks(range(1, num_iter + 1))
        axes[1].set_xticklabels([str(i) for i in range(1, num_iter + 1)])
        axes[1].yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        plt.show()


if __name__ == '__main__':
    s = time.perf_counter()
    sf = Network("SiouxFalls")
    model = ORCS(network=sf,
                 control_intensity=0.1,
                 penalty_param=1,
                 control_potential=1,
                 gap1=1e-5,
                 gap2=1e-5,
                 gap3=1e-5,
                 gap4=1e-5,
                 max_iter1=3000,
                 max_iter2=3000,
                 max_iter3=3000,
                 max_iter4=3000)
    model.run()
    e = time.perf_counter()
    model.plot_basic_data()
    print(f'Total running time: {e-s:.1f}s.')
