import re
import numpy as np
from math import inf


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
        self.cost = (self.fft * (1 + self.alpha * ((self.ue_flow + self.so_flow) / self.capacity) ** self.beta)
         + self.fft * self.alpha * self.beta * ((self.ue_flow + self.so_flow) / self.capacity) ** self.beta)

    def get_cost(self, param) -> float:
        return self.fft * (1 + self.alpha * (param / self.capacity) ** self.beta)

    def get_marginal_cost(self, param) -> float:
        return (self.fft * (1 + self.alpha * (param / self.capacity) ** self.beta)
                + self.fft * self.alpha * self.beta * (param / self.capacity) ** self.beta)


class ODPair:
    def __init__(self, origin, destination, demand):
        self.origin: Node = origin
        self.destination: Node = destination
        self.total_demand: float = demand
        self.ue_demand: float = 0
        self.so_demand: float = 0  # TODO

    def __repr__(self):
        return f'ODPair {self.origin.node_id}->{self.destination.node_id}={self.total_demand}'


class Network:
    def __init__(self, name, sst=1):
        self.name: str = name
        self.sst = sst  # demand sensitivity
        self.Nodes: list[Node] = list()
        self.Links: list[Link] = list()
        self.ODPairs: list[ODPair] = list()
        self.num_node: int = 0
        self.num_link: int = 0
        self.read_network()
        self.read_OD()

    def read_network(self):
        path = f'TransportationNetworks\\{self.name}\\{self.name}_net.txt'
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
            cur_link = Link(index + 1, self.Nodes[int(line[0])], self.Nodes[int(line[1])], float(line[2]),
                        float(line[3]), float(line[4]), float(line[5]), float(line[6]))
            self.Links.append(cur_link)
            self.Nodes[int(line[0])].downstream_link.append(cur_link)
            self.Nodes[int(line[1])].upstream_link.append(cur_link)

    def read_OD(self):
        path = f'TransportationNetworks\\{self.name}\\{self.name}_trips.txt'
        with open(path, 'r', encoding='UTF-8') as file:
            # Process the text file
            lines = file.readlines()
        pattern = re.compile(r'[0-9.]+|Origin')
        data = [pattern.findall(line) for line in lines if len(pattern.findall(line)) != 0]
        total_flow = float(data[1][0])
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
                demand = float(line[2 * i + 1]) * self.sst
                if demand != 0:
                    self.ODPairs.append(ODPair(origin, destination, demand))
        if abs(total_flow - sum([od.total_demand for od in self.ODPairs])) > 1:
            raise ValueError("Demand in the file doesn't match with the total OD flow.")


class MixedEquilibrium:
    def __init__(self, network, ue_demand, so_demand, gap1, gap2, max_iter1, max_iter2):
        self.network: Network = network
        self.ue_demand: list[float] = ue_demand
        self.so_demand: list[float] = so_demand
        self.gap1: float = gap1
        self.gap2: float = gap2
        self.max_iter1: int = max_iter1
        self.max_iter2: int = max_iter2


    def initialize(self):
        for index, od in enumerate(self.network.ODPairs):
            od.ue_demand = self.ue_demand[index]
            od.so_demand = self.so_demand[index]
        for link in self.network.Links[1:]:
            link.update_cost()
            link.update_marginal_cost()

    def main(self):
        pass

    def all_or_nothing_ue(self):
        for link in self.network.Links[1:]:
            link.update_cost()
        for od in self.network.ODPairs:
            origin, destination, ue_demand = od.origin.node_id, od.destination.node_id, od.ue_demand
            shortest_path = self.dijkstra(origin, destination, marginal=False)
            for link in shortest_path:
                link.auxiliary_ue_flow += ue_demand

    def all_or_nothing_so(self):
        for link in self.network.Links[1:]:
            link.update_marginal_cost()
        for od in self.network.ODPairs:
            origin, destination, so_demand = od.origin.node_id, od.destination.node_id, od.so_demand
            shortest_path = self.dijkstra(origin, destination, marginal=False)
            for link in shortest_path:
                link.auxiliary_so_flow += so_demand

    """
    Find the shortest path based on the cost or the marginal cost.
    Set marginal=False for the former or marginal=True for the latter.
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
                return None
            shortest_path.append(temp)
            cur_node = p
        shortest_path.reverse()
        return shortest_path



class ADMM:
    def __init__(self, network, control_intensity, penalty_param, control_potential, gap1, gap2, gap3, gap4, max_iter1,
                 max_iter2, max_iter3, max_iter4):
        self.network = network
        self.control_intensity = control_intensity
        self.penalty_param = penalty_param
        self.control_potential = control_potential
        self.gap1 = gap1
        self.gap2 = gap2
        self.gap3 = gap3
        self.gap4 = gap4
        self.max_iter1 = max_iter1
        self.max_iter2 = max_iter2
        self.max_iter3 = max_iter3
        self.max_iter4 = max_iter4


if __name__ == '__main__':
    sf = Network("SiouxFalls")
    model = ADMM(network=sf,
                 control_intensity=0.1,
                 penalty_param=1,
                 control_potential=1,
                 gap1=1e-4,
                 gap2=1e-4,
                 gap3=1e-4,
                 gap4=1e-4,
                 max_iter1=1000,
                 max_iter2=1000,
                 max_iter3=1000,
                 max_iter4=1000)
    lower = MixedEquilibrium(sf, [100]*528, [100]*528, 1e-4, 1e-4, 1000, 1000)
    lower.initialize()
    print(lower.dijkstra(1, 24, marginal=True))