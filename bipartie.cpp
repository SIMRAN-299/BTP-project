#include <bits/stdc++.h>
using namespace std;

double epsilon = 1e-7;
class Point
{
public:
    double x;
    double y;
    double length;
    string name;
    int id = 0;

    Point(double x, double y, double length, string name)
    {
        this->x = x;
        this->y = y;
        this->length = length;
        this->name = name;
    }
};
class Graph
{
    int V;
    map<int, vector<int>> adj;
    vector<int> color;
    vector<int> visited;

public:
    Graph(int V)
    {
        this->V = V;
        for (int i = 0; i <= V; i++)
        {
            adj[i].clear();
            color.push_back(0);
            visited.push_back(0);
        }
    }
    void addEdge(int u, int v)
    {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void colorNode(int u)
    {
        color[u] = 1;
    }
    void printcolor()
    {
        for (int i = 1; i <= V; i++)
        {
            cout << color[i] << " ";
        }
        cout << endl;
    }
    void dfs(int node, vector<int> &color_1, vector<int> &color_2)
    {

        visited[node] = 1;
        if (color[node] == 0)
        {

            color_1.push_back(node);
        }
        else
        {

            color_2.push_back(node);
        }

        for (int i : adj[node])
        {
            if (!visited[i])
            {
                dfs(i, color_1, color_2);
            }
        }
    }
    vector<int> getNodes()
    {
        vector<int> ans;
        for (int i = 1; i <= V; i++)
        {
            if (!visited[i])
            {
                vector<int> A;
                vector<int> B;
                dfs(i, A, B);

                if (A.size() > B.size())
                {
                    for (int i : A)
                    {
                        ans.push_back(i);
                    }
                }
                else
                {
                    for (int i : B)
                    {
                        ans.push_back(i);
                    }
                }
            }
        }
        return ans;
    }
};
class LabellingProblem
{
public:
    double n;
    double height;
    vector<Point> points;
    LabellingProblem(double n, double height, vector<Point> points)
    {
        this->n = n;
        this->height = height;
        this->points = points;
    }
    bool checkObscureHelper2(pair<double, double> x_range, pair<double, double> y_range, vector<double> point)
    {
        double left = x_range.first;
        double right = x_range.second;
        double down = y_range.first;
        double up = y_range.second;

        double x;
        double y;

        x = point[0];
        y = point[1];
        bool res1 = (left < x || abs(left - x) < epsilon) && (x < right || abs(right - x) < epsilon) && (up > y || abs(up - y) < epsilon) && (y > down || abs(down - y) < epsilon);

        x = point[0];
        y = point[1] + height;
        bool res2 = (left < x || abs(left - x) < epsilon) && (x < right || abs(right - x) < epsilon) && (up > y || abs(up - y) < epsilon) && (y > down || abs(down - y) < epsilon);

        x = point[0] + point[2];
        y = point[1];
        bool res3 = (left < x || abs(left - x) < epsilon) && (x < right || abs(right - x) < epsilon) && (up > y || abs(up - y) < epsilon) && (y > down || abs(down - y) < epsilon);

        x = point[0] + point[2];
        y = point[1] + height;
        bool res4 = (left < x || abs(left - x) < epsilon) && (x < right || abs(right - x) < epsilon) && (up > y || abs(up - y) < epsilon) && (y > down || abs(down - y) < epsilon);

        return res1 || res2 || res3 || res4;
    }
    bool checkObscureHelper1(pair<double, double> x_range, pair<double, double> y_range, vector<double> point)
    {
        double left = x_range.first;
        double right = x_range.second;
        double down = y_range.first;
        double up = y_range.second;

        double x = point[0];
        double y = point[1];

        bool res = (left < x || abs(left - x) < epsilon) && (x < right || abs(right - x) < epsilon) && (up > y || abs(up - y) < epsilon) && (y > down || abs(down - y) < epsilon);
        // bool res=left<x&&x<right&&up>y&&y>down;
        return res;
    }
    bool checkObscure(pair<double, double> x_range, pair<double, double> y_range, double i)
    {
        for (double j = 0; j < n; j++)
        {
            if (i == j)
            {
                continue;
            }
            bool res = checkObscureHelper1(x_range, y_range, {points[j].x, points[j].y});
            if (res == 1)
            {
                return 0;
            }
        }
        return 1;
    }
    vector<vector<double>> processObscuringPoints()
    {

        vector<vector<double>> validLabels;
        for (double i = 0; i < n; i++)
        {
            vector<vector<double>> tempStore;
            double x = points[i].x;
            double y = points[i].y;
            double length = points[i].length;

            pair<double, double> x_range;
            pair<double, double> y_range;
            bool res;

            x_range = {x - length / 2, x + length / 2};
            y_range = {y - height, y};

            res = checkObscure(x_range, y_range, i);
            if (res == 1)
            {

                tempStore.push_back({x_range.first, y_range.first, length, i});
            }

            x_range = {x - length / 2, x + length / 2};
            y_range = {y, y + height};

            res = checkObscure(x_range, y_range, i);
            if (res == 1)
            {
                tempStore.push_back({x_range.first, y_range.first, length, i});
            }
            x_range = {x - length, x};
            y_range = {y - height / 2, y + height / 2};

            res = checkObscure(x_range, y_range, i);
            if (res == 1)
            {
                tempStore.push_back({x_range.first, y_range.first, length, i});
            }
            x_range = {x, x + length};
            y_range = {y - height / 2, y + height / 2};

            res = checkObscure(x_range, y_range, i);
            if (res == 1)
            {
                tempStore.push_back({x_range.first, y_range.first, length, i});
            }
            random_shuffle(tempStore.begin(), tempStore.end());
            for (auto it : tempStore)
            {
                validLabels.push_back(it);
            }
        }

        return validLabels;
    }
    bool static compareH(vector<double> a, vector<double> b)
    {
        return a[0] <= b[0];
    }
    bool static compareV(vector<double> a, vector<double> b)
    {
        return a[1] >= b[1];
    }
    vector<vector<double>> processHorizontalSweep(vector<vector<double>> labels, bool direction)
    {
        sort(labels.begin(), labels.end(), compareH);
        if (direction == 1)
        {
            reverse(labels.begin(), labels.end());
        }
        vector<vector<double>> visited;
        double n = labels.size();
        for (double i = 0; i < n; i++)
        {

            bool res1 = 0, res2 = 0;

            for (auto it : visited)
            {
                pair<double, double> x_range = {labels[i][0], labels[i][0] + labels[i][2]};
                pair<double, double> y_range = {labels[i][1], labels[i][1] + height};

                res1 = checkObscureHelper2(x_range, y_range, it);

                x_range = {it[0], it[0] + it[2]};
                y_range = {it[1], it[1] + height};
                res2 = checkObscureHelper2(x_range, y_range, labels[i]);
                if (res1 || res2)
                {
                    break;
                }
            }
            if (!res1 && !res2)
            {

                visited.push_back(labels[i]);
            }
        }
        return visited;
    }
    vector<vector<double>> processVerticalSweep(vector<vector<double>> labels, bool direction)
    {
        sort(labels.begin(), labels.end(), compareV);
        if (direction == 1)
        {
            reverse(labels.begin(), labels.end());
        }
        vector<vector<double>> visited;
        double n = labels.size();
        for (double i = 0; i < n; i++)
        {
            bool res1 = 0, res2 = 0;
            for (auto it : visited)
            {

                pair<double, double> x_range = {labels[i][0], labels[i][0] + labels[i][2]};
                pair<double, double> y_range = {labels[i][1], labels[i][1] + height};
                res1 = checkObscureHelper2(x_range, y_range, it);

                x_range = {it[0], it[0] + it[2]};
                y_range = {it[1], it[1] + height};
                res2 = checkObscureHelper2(x_range, y_range, labels[i]);
                if (res1 || res2)
                {
                    break;
                }
            }
            if (!res1 && !res2)
            {
                visited.push_back(labels[i]);
            }
        }
        return visited;
    }

    vector<vector<double>> bipartiteMatching(vector<vector<double>> A, vector<vector<double>> B)
    {
        map<vector<double>, int> label_to_number; // mapping each double to an integer
        map<int, vector<double>> number_to_label; // mapping each integer to double

        int size = A.size() + B.size();
        for (auto x : A)
        {
            for (auto y : B)
            {
                if (x == y)
                {
                    size--;
                }
            }
        }
        Graph label_graph(size); // initiallization of graph
        int number = 1;
        for (auto label : A)
        {
            label_to_number[label] = number;
            label_graph.colorNode(number); // coloring nodes of A 1
            number_to_label[number++] = label;
        }
        for (auto label : B)
        {
            if (label_to_number[label] == 0)
            {
                label_to_number[label] = number;

                number_to_label[number++] = label;
            }
        }
        // label_graph.printcolor();
        for (auto labelA : A)
        {
            int node1 = label_to_number[labelA];
            for (auto labelB : B)
            {
                int node2 = label_to_number[labelB];
                pair<double, double> x_range = {labelA[0], labelA[0] + labelA[2]};
                pair<double, double> y_range = {labelA[1], labelA[1] + height};
                bool res1 = checkObscureHelper2(x_range, y_range, labelB);

                x_range = {labelB[0], labelB[0] + labelB[2]};
                y_range = {labelB[1], labelB[1] + height};
                bool res2 = checkObscureHelper2(x_range, y_range, labelA);
                if (res1 || res2)
                {

                    label_graph.addEdge(node1, node2); // add edge between if labels are intersecting
                }
            }
        }
        vector<int> Nodes = label_graph.getNodes();
        vector<vector<double>> ans;
        for (int i : Nodes)
        {
            ans.push_back(number_to_label[i]);
        }

        return ans;
    }
    bool isSafeForIndependentSet(
        int vertex,
        set<int> &tempSolutionSet, map<pair<int, int>, int> &edges)
    {
        for (auto iter : tempSolutionSet)
        {
            if (edges.find(make_pair(iter, vertex)) != edges.end())
            {
                return false;
            }
        }
        return true;
    }
    void ExhaustiveSearch(int currV, int setSize, set<int> &tempSolutionSet, int &ans, vector<int> &vertices, map<pair<int, int>, int> &edges)
    {
        if (currV == setSize)
        {
            int temp = tempSolutionSet.size();
            ans = max(ans, temp);
            return;
        }
        for (int i = currV; i <= setSize; i++)
        {
            if (isSafeForIndependentSet(
                    vertices[i - 1],
                    tempSolutionSet, edges))
            {
                tempSolutionSet
                    .insert(vertices[i - 1]);
                ExhaustiveSearch(
                    i + 1,
                    setSize,
                    tempSolutionSet, ans, vertices, edges);
                tempSolutionSet
                    .erase(vertices[i - 1]);
            }
        }
    }
    vector<vector<double>> solve()
    {
        vector<vector<double>> labels = processObscuringPoints(); // valid labels
        int V = labels.size();
        map<pair<int, int>, int> edges;
        vector<int> vertices;
        int _ans=INT_MIN;
        for (int i = 0; i < V; i++)
        {
            vertices.push_back(i);
        }
        for (int i = 0; i < V; i++)
        {
            for (int j = 0; j < V; j++)
            {
                if (i == j)
                {
                    continue;
                }
                bool res1 = 0, res2 = 0;
                pair<double, double> x_range = {labels[i][0], labels[i][0] + labels[i][2]};
                pair<double, double> y_range = {labels[i][1], labels[i][1] + height};
                res1 = checkObscureHelper2(x_range, y_range, labels[j]);

                x_range = {labels[j][0], labels[j][0] + labels[j][2]};
                y_range = {labels[j][1], labels[j][1] + height};
                res2 = checkObscureHelper2(x_range, y_range, labels[i]);
                if (res1 || res2)
                {
                    edges[make_pair(i, j)] = 1;
                    edges[make_pair(j, i)] = 1;
                }
            }
        }
        set<int> tempSolutionSet;

        ExhaustiveSearch(0, V, tempSolutionSet, _ans, vertices, edges);
        cout << _ans << "----" << "Exhaustive Search" << endl;

        vector<vector<double>> H1 = processHorizontalSweep(labels, 0); // left to right
        vector<vector<double>> H2 = processHorizontalSweep(labels, 1); // right to left
        vector<vector<double>> V1 = processVerticalSweep(labels, 0);   // top to bottom
        vector<vector<double>> V2 = processVerticalSweep(labels, 1);   // bottom to top
        vector<vector<vector<double>>> A = {H1, H2, V1, V2};

        vector<vector<double>> ans = {};
        for (int i = 0; i < 4; i++)
        {
            cout << A[i].size() << endl;
            ans = bipartiteMatching(A[i], ans);
        }
        cout << ans.size() << endl;
        // for(auto i:ans){
        //             cout<<i[0]<<" "<<i[1]<<" "<<points[i[3]].name<<" "<<i[2]<<endl;
        // }

        return ans;
    }
    void verifyOverlapTest(vector<vector<double>> &labels)
    {
        for (int i = 0; i < labels.size(); i++)
        {
            for (int j = 0; j < labels.size(); j++)
            {
                if (i == j)
                {
                    continue;
                }

                bool res1 = 0, res2 = 0;
                pair<double, double> x_range = {labels[i][0], labels[i][0] + labels[i][2]};
                pair<double, double> y_range = {labels[i][1], labels[i][1] + height};
                res1 = checkObscureHelper2(x_range, y_range, labels[j]);

                x_range = {labels[j][0], labels[j][0] + labels[j][2]};
                y_range = {labels[j][1], labels[j][1] + height};
                res2 = checkObscureHelper2(x_range, y_range, labels[i]);
                if (res1 || res2)
                {
                    cout << "failed!" << endl;
                    return;
                }
            }
        }

        cout << "passed!" << endl;
        return;
    }
    void generateSolutionFile(vector<vector<double>> &labels, double height)
    {
        double shift_x = 1e6, shift_y = 1e6;
        for (int i = 0; i < labels.size(); i++)
        {
            shift_x = min(shift_x, labels[i][0]);
            shift_y = min(shift_y, labels[i][1]);
        }
        fstream myFile;
        myFile.open("output.txt", ios::out); // write
        if (myFile.is_open())
        {
            myFile << to_string(height) << "\n";
            for (int i = 0; i < labels.size(); i++)
            {
                double x = labels[i][0] - shift_x;
                double y = labels[i][1] + height - shift_y;
                double len = labels[i][2];
                string line = to_string(x) + " " + to_string(y) + " " + to_string(len) + "\n";
                myFile << line;
            }
            myFile.close();
        }
    }
};

int main()
{
    srand(time(0));
    fstream myFile;
    vector<Point> points;
    myFile.open("input.txt", ios::in); // read
    if (myFile.is_open())
    {
        string line;
        while (getline(myFile, line))
        {
            istringstream ss(line);
            string x_pos, y_pos, name;
            double x, y, len;
            ss >> x_pos;
            ss >> y_pos;
            getline(ss, name);

            x = stod(x_pos);
            y = stod(y_pos);
            name = name.substr(2, name.size() - 3);
            len = name.size();
            Point p(x * 100, y * 100, len, name);
            points.push_back(p);
        }
        myFile.close();
    }
    int n = points.size();
    cout << n << endl;
    int h = 4;
    LabellingProblem obj(n, h, points);
    vector<vector<double>> ans = obj.solve(); // solve
    obj.verifyOverlapTest(ans);               // test overlap
    // obj.generateSolutionFile(ans,h);
    return 0;
}

// int main()
// {
//     srand(time(0));
//     Point p1(3, 2, 1,"A");
//     Point p2(5.5, 4.5, 2,"B");
//     Point p3(7, 1.5, 3.5,"C");
//     Point p4(3, 1.8, 1,"D");
//     Point p5(3, 1.8, 1,"E");
//     Point p6(6, 4.5, 1,"F");
//     Point p7(6, 4.5, 1,"G");
//     Point p8(7.5, 1.5, 1,"H");
//     Point p9(7.5, 1.5, 1,"I");

//     Point p10(20, 21, 1.2,"J");
//     Point p11(23, 21, 1.2,"K");
//     Point p12(21.5, 23, 1,"L");
//     Point p13(21.5, 23.5, 1,"M");
//     Point p14(21.5, 23.5, 1,"N");
//     Point p15(23.5, 21, 1.2,"O");
//     Point p16(23.5, 21, 1.2,"P");
//     Point p17(19.5, 21, 1,"Q");
//     Point p18(19.5, 21, 1,"R");

//     Point p19(49,43,4,"S");
//     Point p20(45.5,44.5, 2,"T");
//     Point p21(47, 41.5, 3.5,"U");
//     Point p22(50, 43, 1,"V");
//     Point p23(50, 43, 1,"W");
//     Point p24(46, 44.5, 1,"X");
//     Point p25(46, 44.5, 1,"Y");
//     Point p26(47.5, 41.5, 1,"Z");
//     Point p27(47.5, 41.5, 1,"a");

//     Point p28(71.5, 70, 1,"b");
//     Point p29(73, 71, 1.2,"c");
//     Point p30(70, 71, 1,"d");
//     Point p31(71.5, 69.5, 1,"f");
//     Point p32(71.5, 69.5, 1,"g");
//     Point p33(73.5, 71, 1.2,"h");
//     Point p34(73.5, 71, 1.2,"l");
//     Point p35(69.5, 71, 1,"m");
//     Point p36(69.5, 71, 1,"n");

//     vector<Point> points = {p1,p2,p3,p4,p5,p6,p7,p8,p9,
//     p10,p11,p12,p13,p14,p15,p16,p17,p18,
//     p19,p20,p21,p22,p23,p24,p25,p26,p27,
//     p28,p29,p30,p31,p32,p33,p34,p35,p36};
//     LabellingProblem obj(36, 2, points);
//     obj.solve();
//     return 0;
// }
