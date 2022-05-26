#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

struct man {
    int index, x, y, z;
};

template<typename Type>
void MergeSort(Type* buf, size_t l, size_t r, bool (*comp)(Type l, Type r));
template<typename Type>
static void merge(Type* buf, size_t l, size_t r, size_t m, bool (*comp)(Type l, Type r));

template<typename Type>
void MergeSort(std::vector<Type>& buf, size_t l, size_t r, bool (*comp)(Type l, Type r))
{

    if (l >= r) return;

    size_t m = (l + r) / 2;


    MergeSort(buf, l, m, comp);
    MergeSort(buf, m + 1, r, comp);
    merge(buf, l, r, m, comp);
}

template<typename Type>
static void merge(std::vector<Type>& buf, size_t l, size_t r, size_t m, bool (*comp)(Type l, Type r))
{
    if (l >= r || m < l || m > r) return;
    if (r == l + 1 && comp(buf[r], buf[l])) { //buf[l] > buf[r]
        std::swap(buf[l], buf[r]);
        return;
    }

    std::vector<Type> tmp(&buf[l], &buf[r] + 1);

    for (size_t i = l, j = 0, k = m - l + 1; i <= r; ++i) {
        if (j > m - l) {
            buf[i] = tmp[k++];
        }
        else if (k > r - l) {
            buf[i] = tmp[j++];
        }
        else {
            buf[i] = (comp(tmp[j], tmp[k])) ? tmp[j++] : tmp[k++]; //tmp[j] < tmp[k]
        }
    }
}

bool Comp_for_x(man l, man r) {
    if (l.x < r.x) return true;
    return false;
}

bool Comp_for_y(man l, man r) {
    if (l.y < r.y) return true;
    return false;
}

bool Comp_for_answer(man l, man r) {
    if (l.index < r.index) return true;
    return false;
}

void build(const std::vector <man>& a, int v, int tl, int tr, std::vector <int>& tree) {
    if (tl == tr)
        tree[v] = a[tl].z;
    else {
        int tm = (tl + tr) >> 1;
        build(a, v * 2, tl, tm, tree);
        build(a, v * 2 + 1, tm + 1, tr, tree);
        tree[v] = std::max(tree[v * 2], tree[v * 2 + 1]);
    }
}

void update(int v, int tl, int tr, int pos, int new_val, std::vector<int>& tree) {
    if (tl == tr)
        tree[v] = new_val;
    else {
        int tm = (tl + tr) >> 1;
        if (pos <= tm)
            update(v * 2, tl, tm, pos, new_val, tree);
        else
            update(v * 2 + 1, tm + 1, tr, pos, new_val, tree);
        tree[v] = std::max(tree[v * 2], tree[v * 2 + 1]);
    }
}

int find_max(int v, int tl, int tr, int l, int r, std::vector <int>& tree) {
    if (l > r) return 0;
    if (l == tl && r == tr) return tree[v];
    int tm = (tl + tr) >> 1;
    return std::max(find_max(v * 2, tl, tm, l, std::min(r, tm), tree),
        find_max(v * 2 + 1, tm + 1, tr, std::max(l, tm + 1), r, tree));
}

int main()
{
    std::ifstream inp("input.txt");
    std::ofstream out("output.txt");

    int n, max = 0;
    inp >> n;

    std::vector<man> mas(n);
    std::vector<man> sort_by_x(n);
    std::vector<man> sort_by_y(n);
    std::vector<int> place_in_y_sort(n);
    std::vector<int> tree(n * 4);
    std::vector<man> answer;
    for (int i = 0; i < n; i++) {
        inp >> mas[i].x >> mas[i].y >> mas[i].z;
        mas[i].index = i + 1;
        sort_by_x[i] = mas[i];
        sort_by_y[i] = mas[i];
    }

    MergeSort(sort_by_x, 0, sort_by_x.size() - 1, Comp_for_x);
    MergeSort(sort_by_y, 0, sort_by_y.size() - 1, Comp_for_y);
    //std::sort(sort_by_x.begin(), sort_by_x.end(), Comp_for_x);
    //std::sort(sort_by_y.begin(), sort_by_y.end(), Comp_for_y);

    build(sort_by_y, 1, 0, sort_by_y.size() - 1, tree);

    for (int i = 0; i < n; i++) {
        place_in_y_sort[sort_by_y[i].index - 1] = i;
    }

    for (int i = 0; i < n; i++) {
        man temp = sort_by_x[i];

        int t = 0;
        while (i + t < n && sort_by_x[i + t].x == sort_by_x[i].x) {
            update(1, 0, n - 1, place_in_y_sort[sort_by_x[i + t].index - 1], 0, tree);
            t++;
        }

        do {
            temp = sort_by_x[i];
            int left = place_in_y_sort[temp.index - 1];
            do {
                left++;
            } while (left < n && sort_by_y[left - 1].y == sort_by_y[left].y);

            if (left < n && sort_by_y[left].y != temp.y) {
                int max_z = find_max(1, 0, n - 1, left, n - 1, tree);
                if (max_z <= temp.z) {
                    answer.push_back(temp);
                }
            }
            else {
                answer.push_back(temp);
            }
            i++;
        } while (i < n && sort_by_x[i - 1].x == sort_by_x[i].x);
        i--;

    }

    MergeSort(answer, 0, answer.size() - 1, Comp_for_answer);
    //std::sort(answer.begin(), answer.end(), Comp_for_answer);

    for (int i = 0; i < answer.size(); i++) {
        out << answer[i].index << ' ';
    }

    return 0;


}
