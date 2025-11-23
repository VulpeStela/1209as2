from collections import defaultdict

class OpenMeanderError(Exception):
    pass

class DyckWordError(Exception):
    pass

class OpenMeander:
    def __init__(self, *args):
        if len(args) < 2:
            raise OpenMeanderError('Not a permutation of 1,..., n for some n >= 2')
        
        sorted_args = sorted(list(args))
        n = len(args)
        if sorted_args != list(range(1, n + 1)):
            raise OpenMeanderError('Not a permutation of 1,..., n for some n >= 2')
        
        self.points = list(args)
        self.n = n
        
        self.upper_edges = []
        self.lower_edges = []
        
        for i in range(n - 1):
            u, v = self.points[i], self.points[i+1]
            edge = tuple(sorted((u, v)))
            if i % 2 == 0:
                self.upper_edges.append(edge)
            else:
                self.lower_edges.append(edge)
                
        if self._check_intersections(self.upper_edges) or self._check_intersections(self.lower_edges):
             raise OpenMeanderError('Does not define an open meander')

        self.extended_dyck_word_for_upper_arches = self._generate_dyck_word(self.upper_edges)
        self.extended_dyck_word_for_lower_arches = self._generate_dyck_word(self.lower_edges)

    def _check_intersections(self, edges):
        count = len(edges)
        for i in range(count):
            for j in range(i + 1, count):
                a, b = edges[i]
                c, d = edges[j]
                if (a < c < b < d) or (c < a < d < b):
                    return True
        return False

    def _generate_dyck_word(self, edges):
        connection_map = {}
        for u, v in edges:
            connection_map[u] = ('start', v)
            connection_map[v] = ('end', u)
        result = []
        for k in range(1, self.n + 1):
            if k not in connection_map:
                result.append('1')
            else:
                role, _ = connection_map[k]
                if role == 'start':
                    result.append('(')
                else:
                    result.append(')')
        return ''.join(result)

    def draw(self, filename, scale=1):
        s = float(scale)
        
        
        
        lines = []
        lines.append(r'\documentclass[10pt]{article}')
        lines.append(r'\usepackage{tikz}')
        lines.append(r'\usepackage[margin=0cm]{geometry}')
        lines.append(r'\pagestyle{empty}')
        lines.append('')
        lines.append(r'\begin{document}')
        lines.append('')
        lines.append(r'\vspace*{\fill}')
        lines.append(r'\begin{center}')
        lines.append(f'\\begin{{tikzpicture}}[x={fmt_float(s)}cm, y={fmt_float(s)}cm, very thick]')

        x_min = 0
        x_max = self.n + 1
        lines.append(f'\\draw ({x_min},0) -- ({x_max},0);')

        upper_connected = {pt for edge in self.upper_edges for pt in edge}
        lower_connected = {pt for edge in self.lower_edges for pt in edge}
        strut_len = 0.5

        for i in range(self.n):
            u = self.points[i]
            if u not in upper_connected:
                lines.append(f'\\draw ({fmt_int(u)},0) -- ({fmt_int(u)},{fmt_float(strut_len)});')
            if u not in lower_connected:
                lines.append(f'\\draw ({fmt_int(u)},0) -- ({fmt_int(u)},-{fmt_float(strut_len)});')

            if i < self.n - 1:
                v = self.points[i+1]
                radius = (v - u) / 2.0
                is_upper = (i % 2 == 0)
                if is_upper == (radius > 0):
                    start_angle = 180
                else:
                    start_angle = -180
                lines.append(f'\\draw ({fmt_int(u)},0) arc[start angle={start_angle}, end angle=0, radius={fmt_float(radius)}];')

        lines.append(r'\end{tikzpicture}')
        lines.append(r'\end{center}')
        lines.append(r'\vspace*{\fill}')
        lines.append('')
        lines.append(r'\end{document}')
        lines.append('')

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))


class DyckWord:
    def __init__(self, word: str):
        if len(word) == 0:
            raise DyckWordError('Expression should not be empty')
        if any(c not in '()' for c in word):
            raise DyckWordError('Expression can only contain \'(\' and \')\'')
            
        self.word = word
        self.arches = [] 
        
        stack = []
        current_depth = 0
        
        for i, char in enumerate(word):
            if char == '(':
                stack.append((i, current_depth))
                current_depth += 1
            elif char == ')':
                if not stack:
                    raise DyckWordError('Unbalanced parentheses in expression')
                start_index, depth = stack.pop()
                self.arches.append({'start': start_index, 'end': i, 'depth': depth})
                current_depth -= 1
                
        if stack:
            raise DyckWordError('Unbalanced parentheses in expression')
            
        self._calculate_properties()

    def _calculate_properties(self):
        sorted_by_depth = sorted(self.arches, key=lambda x: x['depth'], reverse=True)
        children_map = defaultdict(list)
        for arch in self.arches:
            u, v = arch['start'], arch['end']
            potential_children = [a for a in self.arches if a['depth'] == arch['depth'] + 1 and a['start'] > u and a['end'] < v]
            children_map[(u, v)] = sorted(potential_children, key=lambda x: x['start'])

        for arch in sorted_by_depth:
            u, v = arch['start'], arch['end']
            children = children_map[(u, v)]
            span_len = v - u
            
            profile = [0] * span_len
            
            for child in children:
                c_u, c_v = child['start'], child['end']
                c_poly_heights = [h + 1 for h in child['profile']]
                offset = c_u - u
                for k in range(len(c_poly_heights)):
                    profile[offset + k] = c_poly_heights[k]
                    
            # Gap Filling
            prof = list(profile)
            non_zero_indices = [i for i, h in enumerate(profile) if h > 0]
            
            if non_zero_indices:
                for i in range(len(profile)):
                    if profile[i] == 0:
                        left_h = 0
                        left_cands = [idx for idx in non_zero_indices if idx < i]
                        if left_cands:
                            left_h = profile[left_cands[-1]]
                        right_h = 0
                        right_cands = [idx for idx in non_zero_indices if idx > i]
                        if right_cands:
                            right_h = profile[right_cands[0]]
                        prof[i] = max(left_h, right_h)
            
            arch['profile'] = prof
            # depth = height - 1
            arch['height'] = max(prof) + 1 if prof else 1

    def report_on_depths(self):
        depth_counts = defaultdict(int)
        for arch in self.arches:
            d = arch['height'] - 1
            depth_counts[d] += 1
            
        sorted_depths = sorted(depth_counts.keys())
        for d in sorted_depths:
            count = depth_counts[d]
            if count == 1:
                print(f'There is 1 arch of depth {d}.')
            else:
                print(f'There are {count} arches of depth {d}.')

    def _draw(self, filename, scale, coloured=False):
        s = float(scale)
        def fmt(val):
            return '{:.1f}'.format(val)
        
        L = len(self.word)
        lines = []
        lines.append(r'\documentclass[10pt]{article}')
        if coloured:
            lines.append(r'\usepackage[dvipsnames]{xcolour}')
        lines.append(r'\usepackage{tikz}')
        lines.append(r'\usepackage[margin=0cm]{geometry}')
        lines.append(r'\pagestyle{empty}')
        lines.append('')
        lines.append(r'\begin{document}')
        lines.append('')
        lines.append(r'\vspace*{\fill}')
        lines.append(r'\begin{center}')
        lines.append(f'\\begin{{tikzpicture}}[x={fmt(s)}cm, y={fmt(s)}cm, very thick]')
        
        lines.append(f'\\draw (-1,0) -- ({L},0);')

        if coloured:
            draw_order = sorted(self.arches, key=lambda x: (-x['height'], x['start']))
        else:
            draw_order = sorted(self.arches, key=lambda x: x['start'])

        colours = ['Red', 'Orange', 'Goldenrod', 'Yellow', 'LimeGreen', 'Green', 'Cyan', 'SkyBlue', 'Blue', 'Purple']
        last_depth = None

        for arch in draw_order:
            u, v = arch['start'], arch['end']
            profile = arch['profile']
            points = []
            points.append(f'({u},0)')
            
            curr_y = 0
            for i in range(len(profile)):
                h = profile[i] + 1
                x = u + i
                if h != curr_y:
                    points.append(f'({x},{h})')
                    curr_y = h
                next_h = -1
                if i + 1 < len(profile):
                    next_h = profile[i+1] + 1
                if h != next_h:
                    points.append(f'({x+1},{h})')
            
            points.append(f'({v},0)')
            path_str = ' -- '.join(points)
            
            depth = arch['height'] - 1
            
            if coloured:
                if depth != last_depth:
                    lines.append(f'% Arches of depth {depth}')
                    last_depth = depth
                
                c_idx = depth % 10
                colour = colours[c_idx]
                lines.append(f'    \\draw[fill={colour}] {path_str};')
            else:
                lines.append(f'% Arch {draw_order.index(arch) + 1}')
                lines.append(f'    \\draw {path_str};')

        lines.append(r'\end{tikzpicture}')
        lines.append(r'\end{center}')
        lines.append(r'\vspace*{\fill}')
        lines.append('')
        lines.append(r'\end{document}')
        lines.append('')

        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    def draw_arches(self, filename, scale=1):
        self._draw(filename, scale, coloured=False)

    def colour_arches(self, filename, scale=1):
        self._draw(filename, scale, coloured=True)

def fmt_float(f):
    return '{:.1f}'.format(f)

def fmt_int(i):
    return str(int(i))

if __name__ == '__main__':
    w = DyckWord('(((((((((((((())))))))))))))')
    w.report_on_depths()
    w.draw_arches('drawn_dyck_word_1.tex', 0.5)
    w.colour_arches('coloured_dyck_word_1.tex', 0.5)
    print('*' * 30)
    w = DyckWord('(()(()(())))')
    w.report_on_depths()
    w.draw_arches('drawn_dyck_word_2.tex')
    w.colour_arches('coloured_dyck_word_2.tex')
    print('*' * 30)
    w = DyckWord('((()())(()(()())))')
    w.report_on_depths()
    w.draw_arches('drawn_dyck_word_3.tex', 0.6)
    w.colour_arches('coloured_dyck_word_3.tex', 0.6)
    print('*' * 30)
    w = DyckWord('((()(()())(()(()(())))((()()))()(()())))')
    w.report_on_depths()
    w.draw_arches('drawn_dyck_word_4.tex', 0.3)
    w.colour_arches('coloured_dyck_word_4.tex', 0.3)
