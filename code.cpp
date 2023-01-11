#include <bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
typedef long long ll;
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
    int exhaustiveSearch(vector<int> component, vector<vector<int>> adj)
    {
        int ans = 0;
        int n = component.size();
        
        for (ll mask = 0; mask < ((1ll) << (n)); mask++)
        {
            vector<int> nodes;
            // cout<<mask<<endl;
            for (ll node = 0; node <= (n - 1); node++)
            {
                if ((1<<node)&(mask))
                {
                    nodes.push_back(node);
                }
            }
            int s = nodes.size();
            
            bool f = 0;
            for (int i = 0; i < s; i++)
            {
                for (int j = 0; j < s; j++)
                {
                    if (adj[component[nodes[i]]][component[nodes[j]]])
                    {
                        f = 1;
                        break;
                    }
                }
                if (f)
                {
                    break;
                }
            }
            if (!f)
            {
                ans = max(ans, s);
            }
        }
        return ans;
    }
    void findComponent(int node, int n, vector<vector<int>> adj, vector<int> &component, vector<int> &vis)
    {
        component.push_back(node);
        vis[node] = 1;
        for (int i = 0; i < n; i++)
        {
            if (adj[node][i] && !vis[i])
            {
                findComponent(i, n, adj, component, vis);
            }
        }
    }
    vector<vector<double>> ExhaustiveFiltering(int &ansExhaustive)
    {
        vector<vector<vector<double>>> tempLabels;
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
            tempLabels.push_back(tempStore);
        }
        
        for(int i=0;i<tempLabels.size();i++){
            bool flag1=1;
            for(int p1=0;p1<tempLabels[i].size();p1++){
                bool flag=1;
                int j=0;
                for(;j<tempLabels.size();j++){
                    if(i==j){
                        continue;
                    }
                    for(int p2=0;p2<tempLabels[j].size();p2++){
                        bool res1=0,res2=0;
                        pair<double, double> x_range = {tempLabels[i][p1][0], tempLabels[i][p1][1] + tempLabels[i][p1][2]};
                        pair<double, double> y_range = {tempLabels[i][p1][1], tempLabels[i][p1][1] + height};

                        res1 = checkObscureHelper2(x_range, y_range, tempLabels[j][p2]);

                        x_range = {tempLabels[j][p2][0], tempLabels[j][p2][0] + tempLabels[j][p2][2]};
                        y_range = {tempLabels[j][p2][1], tempLabels[j][p2][1] + height};
                        res2 = checkObscureHelper2(x_range, y_range, tempLabels[i][p1]);
                        if(res1||res2){
                            flag=0;
                            break;
                        }
                    }
                    if(flag==0){
                        break;
                    }
                }
                if(j==tempLabels.size()){
                    flag1=0;
                    ansExhaustive++;
                    break;
                }
            }   
            if(flag1==1){
                for (auto it : tempLabels[i])
                {
                    validLabels.push_back(it);
                }
            }    
        }
        // cout<<exhaustiveAns<<"---> Exhaustive Ans"<<endl;
        return validLabels;
    }
    vector<vector<double>> solve()
    {
        vector<vector<double>> labels = processObscuringPoints(); // valid labels
        int ansExhaustive = 0;
        vector<vector<double>>elabels=ExhaustiveFiltering(ansExhaustive);
        int V = elabels.size();
        vector<vector<int>> adj(V, vector<int>(V, 0));
        vector<int> vis(V, 0);
        // map<pair<int, int>, int> edges;
        vector<int> vertices;
        int _ans = INT_MIN;
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
                pair<double, double> x_range = {elabels[i][0], elabels[i][0] + elabels[i][2]};
                pair<double, double> y_range = {elabels[i][1], elabels[i][1] + height};
                res1 = checkObscureHelper2(x_range, y_range, elabels[j]);

                x_range = {elabels[j][0], elabels[j][0] + elabels[j][2]};
                y_range = {elabels[j][1], elabels[j][1] + height};
                res2 = checkObscureHelper2(x_range, y_range, elabels[i]);
                if (res1 || res2)
                {
                    adj[i][j] = 1;
                    adj[j][i] = 1;
                    
                }
            }
        }
        
        
        for (int i = 0; i < V; i++)
        {
            if (!vis[i])
            {
                vector<int> component;
                findComponent(i, V, adj, component, vis);
                int x=component.size();
                 cout<<x<<endl;
                if(x>20){
                    int z=exhaustiveSearch(component,adj);
                    cout<<"for 29                     "<< z<<endl;
                // ansExhaustive += exhaustiveSearch(component, adj);
                }
            }
        }
      
        cout << ansExhaustive << "----"
             << "Exhaustive Search" << endl;

        vector<vector<double>> H1 = processHorizontalSweep(labels, 0); // left to right
        vector<vector<double>> H2 = processHorizontalSweep(labels, 1); // right to left
        vector<vector<double>> V1 = processVerticalSweep(labels, 0);   // top to bottom
        vector<vector<double>> V2 = processVerticalSweep(labels, 1);   // bottom to top
        vector<vector<vector<double>>> A = {H1, H2, V1, V2};

        vector<vector<double>> ans = {};
        for (int i = 0; i < 4; i++)
        {

            ans = bipartiteMatching(A[i], ans);
        }
        cout << "Middle point Labelling Heuristic " << ans.size() << endl;

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
    void generateSolutionFile(vector<vector<double>> &labels, double height, string name)
    {
        double shift_x = 1e6, shift_y = 1e6;
        for (int i = 0; i < labels.size(); i++)
        {
            shift_x = min(shift_x, labels[i][0]);
            shift_y = min(shift_y, labels[i][1]);
        }

        fstream myFile;
        string file = "output/";
        int ptr = 0;
        while (ptr < name.length())
        {
            if (name[ptr] == '.')
            {
                name = name.substr(0, ptr);
                break;
            }
            ptr++;
        }
        file += name;
        file += "_output";
        file += ".txt";
        myFile.open(file, ios::out); // write
        if (myFile.is_open())
        {
            myFile << to_string(height) << "\n";
            for (int i = 0; i < labels.size(); i++)
            {
                double x = labels[i][0] /*-shift_x*/;
                double y = labels[i][1] /*+height/*-shift_y*/;
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
    ifstream fin;
    string dir, filepath;
    DIR *dp;
    int num;
    struct dirent *dirp;
    struct stat filestat;

    // getline( cin, dir );  // gets everything the user ENTERs
    dir = "testCases";
    dp = opendir(dir.c_str());
    if (dp == NULL)
    {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }
    int count = 0;

    while ((dirp = readdir(dp)))
    {

        filepath = dir + "/" + dirp->d_name;
        string inputFileName = dirp->d_name;
        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat(filepath.c_str(), &filestat))
            continue;
        if (S_ISDIR(filestat.st_mode))
            continue;

        fstream myFile;
        vector<Point> points;

        myFile.open(filepath, ios::in); // read
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
                Point p(x, y, len, name);
                points.push_back(p);
            }
            myFile.close();
        }
        int n = points.size();
        cout << n << endl;
        double h = 4;
        LabellingProblem obj(n, h, points);
        vector<vector<double>> ans = obj.solve(); // solve
        obj.verifyOverlapTest(ans);               // test overlap
        obj.generateSolutionFile(ans, h, inputFileName);
        count++;
    }

    closedir(dp);

    return 0;
}