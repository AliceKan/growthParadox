import numpy as np
import random
import matplotlib.pyplot as plt
from collections import deque
import math

class Cell:
    def __init__(self, infinite = False):
        self.age = 0  # 用于控制细胞寿命
        self.alive = True  # 细胞状态
        self.mature_age=24
        self.infinite=infinite
        self.pmax=float('inf') if infinite else 10
        self.apoptosis_chance=0 if infinite else 0.1


    def age_up(self):
        self.age += 1.6
        if not self.infinite and self.pmax <=0:
            self.alive=False

    def is_mature(self):
        return self.age >= self.mature_age

    def can_divide(self):
        return self.alive and (self.infinite or self.pmax > 0)

class ImmuneCell:
    def __init__(self):
        self.alive = True  # 免疫细胞状态

    def increase_apoptosis(self, cell):
        #增加相邻非免疫细胞的凋亡概率
        if isinstance(cell, Cell):
            if cell.infinite:
                cell.apoptosis_chance += 0.005
            else:
                cell.apoptosis_chance += 0.01

class Lattice:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.grid = [[None for _ in range(width)] for _ in range(height)]
        self.infi_counts = 1
        self.latticehour = 0
       
    def add_cell(self, x, y, infinite=False):
        if self.grid[y][x] is None:
            self.grid[y][x] = Cell(infinite=infinite)

    def add_immune_cell(self, x, y):
        #添加免疫细胞
        if self.grid[y][x] is None and self.is_on_surface(x,y):
            self.grid[y][x] = ImmuneCell()
        else:
             # 位置不在表层，在外部随机空位生成免疫细胞
            random_position = self.find_random_empty_position_outside_cluster()
            if random_position:
                rx, ry = random_position
                self.grid[ry][rx] = ImmuneCell()

    def find_random_empty_position_outside_cluster(self):
        #在细胞团外部找到一个随机空位
        empty_positions = []

        # 找到所有外部的空位
        for y in range(self.height):
            for x in range(self.width):
                if self.grid[y][x] is None and self.is_connected_to_boundary(x, y):
                    empty_positions.append((x, y))

        if empty_positions:
            # 随机选择一个空位
            return random.choice(empty_positions)
        else:
            return None

    def is_on_surface(self, x, y):
        #判断细胞是否在表层（至少有一个相邻的空格是外部空格）
        neighbors = self.get_neighbors(x, y)
        for nx, ny in neighbors:
            if self.grid[ny][nx] is None and self.is_connected_to_boundary(nx, ny):
                return True
        return False

    def is_connected_to_boundary(self, x, y, visited=None):
        #检查空格是否连通到边界，用于判断是否是表层
        if visited is None:
            visited = set()

        # 如果当前格子已经访问过，直接返回False，避免死循环
        if (x, y) in visited:
            return False
        visited.add((x, y))

        # 如果空格位于边界，说明它连接到外部
        if x == 0 or y == 0 or x == self.width - 1 or y == self.height - 1:
            return True

        # 递归检查相邻的空格
        neighbors = self.get_neighbors(x, y)
        for nx, ny in neighbors:
            if self.grid[ny][nx] is None:  # 如果相邻格子是空的
                if self.is_connected_to_boundary(nx, ny, visited):
                    return True
        return False


    def has_nearby_non_immune_cells(self, x, y):
        #检查免疫细胞周围是否有非免疫细胞
        neighbors = self.get_neighbors(x, y)
        for nx, ny in neighbors:
            cell = self.grid[ny][nx]
            if isinstance(cell, Cell) and not isinstance(cell, ImmuneCell):
                return True,ny,nx
        return False,y,x

    def find_nearest_non_immune_cell(self, x, y):
        #找到距离 (x, y) 最近的非免疫细胞
        nearest_cell = None
        nearest_distance = float('inf')

        for ny in range(self.height):
            for nx in range(self.width):
                cell = self.grid[ny][nx]
                if isinstance(cell, Cell) and not isinstance(cell, ImmuneCell):
                    distance = math.sqrt((nx - x) ** 2 + (ny - y) ** 2)
                    if distance < nearest_distance:
                        nearest_distance = distance
                        nearest_cell = (nx, ny)
        
        return nearest_cell

    def step(self, ifmove):
        queue=deque()
        # 遍历 lattice
        for y in range(self.height):
            for x in range(self.width):
                if self.grid[y][x]:
                    queue.append((x,y))
        queue_list=list(queue)
        random.shuffle(queue_list)
        queue=deque(queue_list)

        while queue:
            x,y=queue.popleft()
            cell = self.grid[y][x]
            if isinstance(cell, Cell):
                cell.age_up()
                #print(f"Cell at ({x}, {y}) age: {cell.age} hours, max divisions: {cell.pmax}")
                if cell.apoptosis_chance > 1:
                    cell.alive = False
                if not cell.alive:
                    self.grid[y][x] = None  # 细胞死亡
                    continue
                    
                else:
                    move_result=self.move_cell(x, y)# 尝试移动细胞
                    if move_result is not None:
                        if ifmove:
                            new_x, new_y = move_result
                        else:
                            new_x, new_y = x,y
                        if self.grid[new_y][new_x] is None:
                            self.grid[new_y][new_x] = cell # 移动细胞
                            self.grid[y][x] = None  # 原位置清空
                            #print(f"cell moved to {new_x} {new_y}")

                        if self.grid[new_y][new_x] and self.grid[new_y][new_x].is_mature() and self.grid[new_y][new_x].can_divide():
                            if not self.grid[new_y][new_x].infinite and random.random()<self.grid[new_y][new_x].apoptosis_chance:
                                #self.grid[new_y][new_x].alive=False
                                self.grid[new_y][new_x]=None
                                if self.latticehour > 400:
                                    self.add_immune_cell(x, y)
                                    #print(f"celldeathandimmune")
                                continue
                            else:
                                #print(f"Replicating cell at ({new_x}, {new_y})")
                                self.replicate_cell(new_x, new_y)  # 尝试在新位置复制细胞
            elif isinstance(cell, ImmuneCell):
                iftarget,target_y,target_x=self.has_nearby_non_immune_cells(x, y)
                if not iftarget:
                    nearest_non_immune = self.find_nearest_non_immune_cell(x, y)
                    if nearest_non_immune is not None:
                        nx, ny = nearest_non_immune

                        # 移动方向
                        dx = nx - x
                        dy = ny - y

                        # 归一化方向到单位步长
                        if dx != 0:
                            dx = int(dx / abs(dx))
                        if dy != 0:
                            dy = int(dy / abs(dy))

                        new_x = x + dx
                        new_y = y + dy

                        # 如果新位置为空，可以移动
                        if  self.grid[new_y][new_x] is None:
                            self.grid[new_y][new_x] = self.grid[y][x]  # 移动细胞
                            self.grid[y][x] = None  # 清空原位置
                        
                else:
                    target_cell = self.grid[target_y][target_x]
                    if self.latticehour > 400:
                        cell.increase_apoptosis(target_cell)
                        #print(target_cell.apoptosis_chance)
 
        self.latticehour += 1.6


    
    


    def move_cell(self, x, y):
        # 定义可能的移动方向
        directions = [
            (0, -1), (0, 1), (-1, 0), (1, 0),  # 上下左右
            (-1, -1), (-1, 1), (1, -1), (1, 1)  # 四个对角线方向
        ]
        random.shuffle(directions)
        
        for dx, dy in directions:
            new_x = x + dx
            new_y = y + dy
            if 0 <= new_x < self.width and 0 <= new_y < self.height and self.grid[new_y][new_x] is None:
                return new_x, new_y  # 如果方向有效，移动到新的位置
        return None  # 如果没有有效方向，保持在原位置

    def replicate_cell(self, x, y):
        #if self.grid[y][x] is None:
            #print(f"Error: Attempted to replicate from an empty cell at ({x}, {y})")
            #return
        # 随机尝试复制细胞到相邻空位置
        directions = [
            (0, -1), (0, 1), (-1, 0), (1, 0),  # 上下左右
            (-1, -1), (-1, 1), (1, -1), (1, 1)  # 四个对角线方向
        ]
        random.shuffle(directions)
        parent_cell = self.grid[y][x]
        for dx, dy in directions:
            new_x, new_y = x + dx, y + dy
            if 0 <= new_x < self.width and 0 <= new_y < self.height:
                if self.grid[new_y][new_x] is None:
                    if parent_cell.infinite:
                        new_infinite=random.random()<0.01
                        new_cell=Cell(infinite=new_infinite)
                        if new_infinite:
                            self.infi_counts += 1
                            #print(f"infinitecell+1")
                    else:
                        new_cell = Cell()
                        new_cell.pmax = parent_cell.pmax - 1  # 分裂上限减1
                        parent_cell.pmax -= 1  # 父细胞分裂上限减1
                        
                    self.grid[new_y][new_x] = new_cell
                    parent_cell.age = 0  # 父细胞重新开始计时，需要24小时成熟
                    break  # 复制一次后结束

    def get_neighbors(self, x, y):
        #获取相邻的8个邻居格子
        directions = [
            (0, -1), (0, 1), (-1, 0), (1, 0),
            (-1, -1), (-1, 1), (1, -1), (1, 1)
        ]
        neighbors = [
            (x + dx, y + dy) for dx, dy in directions
            if 0 <= x + dx < self.width and 0 <= y + dy < self.height
        ]
        return neighbors


    def get_grid_state(self):
        # 返回一个二维数组，表示 lattice 当前状态
        return np.array([[5 if isinstance(cell, ImmuneCell) else 4 if cell and cell.infinite  else 3 if cell and cell.pmax > 10 else 2 if cell and cell.pmax > 5 else 1 if cell else 0 for cell in row] for row in self.grid])

    def get_cell_count(self):
        return sum(cell is not None and isinstance(cell,Cell) for row in self.grid for cell in row)

    def get_non_immune_cell_count(self):
        #计算当前非免疫细胞的数量
        noncount = 0
        for row in self.grid:
            for cell in row:
                if isinstance(cell, Cell) and not isinstance(cell, ImmuneCell):
                    noncount += 1
        return noncount

    
total_hours = 0
for simulation_count in range(20):

    # 初始化 lattice
    lattice = Lattice(100, 100)
    lattice.add_cell(50, 50, infinite=True)

    cell_counts=[]
    nonimmune_counts=[]


    #total_steps=481
    step = 0
    

    # 运行模拟并绘制图像
    #for step in range(total_steps):  # 模拟  步
    while lattice.get_non_immune_cell_count() > 0:
        ifmove=((step+1)%1==0)
        lattice.step(ifmove=ifmove)
        cell_counts.append(lattice.get_cell_count())
        nonimmune_counts.append(lattice.get_non_immune_cell_count())
    
        #-if (step+1)%5==0:

            # 获取当前 lattice 的状态
            #-grid_state = lattice.get_grid_state()
        

            #-cmap = plt.colormaps['coolwarm']
            ##cmap.set_under(color='white')  # 设置背景色为白色
            #-plt.figure(figsize=(6, 6))
            #-plt.imshow(grid_state, cmap=cmap, interpolation='nearest', vmin=0, vmax=5)
            #-plt.title(f'Time: {1.6 * (step + 1):.1f} Hours')
            #-plt.xlabel('X-axis')
            #-plt.ylabel('Y-axis')
            #-plt.colorbar(ticks=[0, 1, 2, 3, 4, 5], label='Cell Type')
            #-plt.show()

        step += 1
    

    total_steps = step
    # 绘制细胞总数随时间变化的图
    #-hours = [4.8 * (i + 1) for i in range(total_steps)]
    #-plt.figure(figsize=(8, 6))
    #-plt.plot(hours, cell_counts, marker='o')
    #-plt.title('Total Cell Count Over Time')
    #-plt.xlabel('Time (hours)')
    #-plt.ylabel('Total Cell Count')
    #-plt.grid(True)
    #-plt.show()
    print("hours:" + str(total_steps*1.6))
    

    #infinite_counts=lattice.infi_counts
    #print(infinite_counts)

    total_hours += total_steps*1.6
avg_hours = total_hours/20
print("avgHour:" + str(avg_hours))
