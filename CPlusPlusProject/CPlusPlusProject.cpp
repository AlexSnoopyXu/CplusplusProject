// CPlusPlusProject.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cMath>
#include <Windows.h>
#include <string>
#include <stack>
#include <queue>
#include <tuple>

using namespace std;

class Node {
public:
    int val;
    vector<Node*> neighbors;

    Node() {
        val = 0;
        neighbors = vector<Node*>();
    }

    Node(int _val) {
        val = _val;
        neighbors = vector<Node*>();
    }

    Node(int _val, vector<Node*> _neighbors) {
        val = _val;
        neighbors = _neighbors;
    }
};

struct ListNode {
	int val;
	ListNode* next;
    ListNode() : val(0), next(nullptr) {}
	ListNode(int x) : val(x), next(NULL) {}
    ListNode(int x, ListNode* next) : val(x), next(next) {}
};

struct TreeNode {
	int val;
	TreeNode* left;
	TreeNode* right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode* left, TreeNode* right) : val(x), left(left), right(right) {}
};

//ListNode* addTwoNumbers(ListNode* l1, ListNode* l2) {
//    ListNode* result;
//    
//    int sum = l1->val + l2->val;
//    ListNode r = sum % 10;
//    int temp = sum / 10;
//    result = &r;
//
//
//    if (l1->next != NULL) {
//        result->next = addTwoNumbers(l1->next, l2->next);
//        result->next->val += temp;
//    }
//    return result;
//}

//int lengthOfLongestSubstring(string s) {
//    map<char, int> charIndexMap;
//    map<char, int> charIndexMap2;
//    unordered_set<char> cSet;
//    int max = 0;
//    int end = -1;
//    int n = s.size();
//    for (int i = 0; i < n; ++i) 
//    {
//        if (i != 0) {
//            // 左指针向右移动一格，移除一个字符
//            cSet.erase(s[i - 1]);
//        }
//        while (end + 1 < n && !cSet.count(s[end + 1])) {
//            // 不断地移动右指针
//            cSet.insert(s[end + 1]);
//            ++end;
//        }
//        // 第 i 到 rk 个字符是一个极长的无重复字符子串
//        max = std::_Max_value(max, end - i + 1);
//    }
//
//    //for (int i = 0; i < s.size(); ++i)
//    //{
//    //    auto c = charIndexMap.find(s[i]);
//    //    // 存在元素，并找到相同位置
//    //    if (c != charIndexMap.end())
//    //    {
//    //        if (charIndexMap.size() > max) 
//    //        {
//    //            max = charIndexMap.size();
//    //        }
//    //        
//    //        i = c->second;
//    //        charIndexMap.clear();
//    //    }
//    //    else
//    //    {
//    //        // 没找到就插入一个
//    //        charIndexMap[s[i]] = i;
//    //    }
//    //}
//
//    //for (int i = 0; i < s.size(); ++i)
//    //{
//    //    auto c = charIndexMap.find(s[i]);
//    //    if (c != charIndexMap.end())
//    //    {
//    //        if (charIndexMap.size() > max)
//    //        {
//    //            max = charIndexMap.size();
//    //        }
//
//    //        // 删除相同项之前的内容
//    //        for (auto cc : charIndexMap)
//    //        {
//    //            if (cc.second > c->second)
//    //            {
//    //                charIndexMap2[cc.first] = cc.second;
//    //            }
//    //        }
//
//    //        charIndexMap2[s[i]] = i;
//    //        charIndexMap.clear();
//    //        charIndexMap.swap(charIndexMap2);
//    //    }
//    //    else
//    //    {
//    //        charIndexMap[s[i]] = i;
//    //    }
//    //}
//
//    //if (charIndexMap.size() > max)
//    //{
//    //    max = charIndexMap.size();
//    //}
//
//    return max;
//}

//string longestPalindrome(string s) {
//    if (s.size() < 2)
//        return s;
//
//    const int n = s.size();
//    int maxLength = 1;
//    int begin = 0;
//    //vector<vector<bool>>  dp(n, vector<bool>(n));
//    vector<bool> dp(((n * n - n) / 2), true);
//
//    //for (int i = 0; i < n; ++i)
//    //{
//    //    dp[i][i] = true;
//    //}
//
//    for (int j = 1; j < n; ++j)
//    {
//        for (int i = 0; i < j; ++i)
//        {
//            /*if (s[i] != s[j]) 
//            {
//                dp[i + j] = false;
//            }
//            else
//            {
//                if (j - i < 3)
//                {
//                    dp[i + j] = true;
//                }
//                else
//                {
//                    dp[i + j] = dp[i + 1 + j - 1];
//                }
//            }*/
//
//            dp[i + j] = (s[i] != s[j]) && (j - i < 3 || dp[i + 1 + j - 1]);
//
//            // 因为j,i最终会相等，即不是一个范围，所以j-i必然为0，j-i+1表示最小长度
//            if (dp[i + j] && (j - i + 1 > maxLength))
//            {
//                maxLength = j - i + 1;
//                begin = i;
//            }
//        }
//    }
//
//    return s.substr(begin, maxLength);
//}

//int rob(vector<int>& nums) {
//    int n = nums.size();
//
//	if (n == 0)
//		return 0;
//    if (n == 1)
//        return nums[0];
//    vector<int> dp(n, 0);
//    dp[0] = nums[0];
//    dp[1] = std::_Max_value(nums[0],nums[1]);
//    for (int i = 2; i < n; ++i)
//    {
//        // 偷当前的房子：前n-2个房子的最大值加当前的钱
//        // 不偷当前房子：前n-1个房子的最大值
//        dp[i] = std::_Max_value(dp[i - 1], (nums[i] + dp[i - 2]));
//    }
//
//    return dp[n - 1];
//}

//int coinChange(vector<int>& coins, int amount) {
//
//    if (coins.size() == 0 || amount == 0)
//    {
//        return -1;
//    }
//
//    vector<int> dp(amount + 1, amount + 1);
//    dp[0] = 0;
//    int n = dp.size();
//    for (int i = 1; i < n; ++i)
//    {
//        // 只增加当前种类一枚币的最小结果
//        for (auto c : coins)
//        {
//            if ((i - c) < 0)
//            {
//                continue;
//            }
//            dp[i] = std::_Min_value(dp[i], 1 + dp[i - c]);
//        }
//    }
//
//    return (dp[amount] == amount + 1) ? -1 : dp[amount];
//}

//int numSquares(int n) {
//    vector<int> dp(n + 1, n + 1);
//    int length = sqrt(n);
//    for (int i = 0; i <= length; ++i)
//    {
//        dp[i * i] = 1;
//    }
//
//    if (dp[n] == 1) {
//        return 1;
//    }
//
//    for (int i = 0; i <= n; ++i)
//    {
//        int temp = sqrt(i);
//        for (int j = 1; j <= temp; ++ j) 
//        {
//            dp[i] = std::_Min_value(dp[i], dp[i - j * j] + 1);
//        }
//    }
//
//    return dp[n];
//}

//int integerBreak(int n) {
//    if (n == 2)
//    {
//        return 1;
//    }
//
//    vector<int> dp(n + 1, 1);
//    for (int i = 3; i <= n; ++i)
//    {
//        for(int j = 1;j <i;++j)
//        {
//            dp[i] = max(dp[i], max(j * (i - j), j * dp[i - j]));
//        }
//    }
//
//    return dp[n];
//}

//bool canPartition(vector<int>& nums) {
//    int n = nums.size();
//    if (n < 2) 
//    {
//        return false;
//    }
//
//    if (n == 2)
//    {
//        return nums[0] == nums[1];
//    }
//
//    int half = n / 2;
//    vector<int> sum(half * n, 0);
//    for (int k = 0; k < half; ++k) 
//    {
//        for (int m = 0; m < n; ++m)
//        {
//            int left = nums[m];
//            for (int i = 1; i <= k; ++i)
//            {
//                left += nums[i];
//            }
//            int right = 0;
//            for (int i = k + 1; i < n; ++i)
//            {
//                right += nums[i];
//            }
//
//            if (left == right)
//            {
//                return true;
//            }
//        }
//        
//    }
//    
//    return false;
//}

//string convert(string s, int numRows) {
//    int n = s.size();
//    if (n <= 1 || n <= numRows || numRows == 1)
//    {
//        return s;
//    }
//
//    string result;
//    int base = numRows + (numRows - 2);
//    for (int i = 0; i < numRows; ++i)
//    {
//        bool b = (i != 0 && i != (numRows - 1));
//        result.push_back(s[i]);
//
//        int j = i + base - (2 * i);
//        while (j < n)
//        {
//            if (j < 0)
//            {
//                break;
//            }
//
//            if (b)
//            {
//                result.push_back(s[j]);
//            }
//
//            int k = j + (2 * i);
//            if (k < 0 || k >= n)
//            {
//                break;
//            }
//            result.push_back(s[k]);
//
//            j += base;
//        }
//    }
//
//    return result;
//}

//int reverse(int x) {
//    
//    if (x == 0)
//    {
//        return x;
//    }
//
//    vector<int> result;
//
//    bool negative = x < 0;
//    int n = abs(x);
//    while (n > 0)
//    {
//        result.push_back((n % 10));
//
//        n /= 10;
//    }
//
//    uint32_t num = 0;
//
//    for (int i = result.size() - 1; i >= 0 ; --i)
//    {
//        if (i > 9 || (i == 9 && result[(result.size() - 1 - i)] > 2))
//        {
//            return 0;
//        }
//        num += pow(10 ,i) * result[(result.size() - 1 - i)];
//    }
//
//    int m = (int)num;
//
//    if (m < 0)
//    {
//        return 0;
//    }
//
//    m *= negative ? -1 : 1;
//
//    return m;
//}

//int minimumTotal(vector<vector<int>>& triangle) {
//    int n = triangle.size();
//
//    vector<vector<int>> dp(n);
//    dp[0].push_back(triangle[0][0]);
//    for (int i = 1; i < n; ++i)
//    {
//        int m = triangle[i].size();
//        for (int j = 0; j < m; ++j)
//        {
//            int left = triangle[i][j];
//            int right = triangle[i][j];
//            if (j - 1 >= 0) 
//            {
//                left += dp[i - 1][j - 1];
//            }
//            else
//            {
//                left = 0x0f0f;
//            }
//
//            if (j < m - 1)
//            {
//                right += dp[i - 1][j];
//            }
//            else
//            {
//                right = 0x0f0f;
//            }
//
//            if (left == 0x0f0f)
//            {
//                dp[i].push_back(right);
//            }
//            else if (right == 0x0f0f)
//            {
//                dp[i].push_back(left);
//            }
//            else
//            {
//                dp[i].push_back(min(left, right));
//            } 
//        }
//    }
//
//    int k = dp[dp.size() - 1].size();
//    int minnum = dp[dp.size() - 1][0];
//    for (int i = 0; i < k; ++i)
//    {
//        if (dp[dp.size() - 1][i] < minnum)
//        {
//            minnum = dp[dp.size() - 1][i];
//        }
//    }
//
//    return minnum;
//}

//int minimumTotal(vector<vector<int>>& triangle) {
//    int n = triangle.size();
//
//    if (n == 0)
//    {
//        return 0;
//    }
//
//    if (n == 1)
//    {
//        return triangle[0][0];
//    }
//
//    vector<vector<int>> dp(n);
//    dp[0].push_back(triangle[0][0]);
//    int minnum = 0x0fff;
//    for (int i = 1; i < n; ++i)
//    {
//        int m = triangle[i].size();
//        for (int j = 0; j < m; ++j)
//        {
//            int left = triangle[i][j];
//            int right = triangle[i][j];
//            if (j - 1 >= 0) 
//            {
//                left += dp[i - 1][j - 1];
//            }
//            else
//            {
//                left = 0x0f0f;
//            }
//
//            if (j < m - 1)
//            {
//                right += dp[i - 1][j];
//            }
//            else
//            {
//                right = 0x0f0f;
//            }
//
//            if (left == 0x0f0f)
//            {
//                dp[i].push_back(right);
//            }
//            else if (right == 0x0f0f)
//            {
//                dp[i].push_back(left);
//            }
//            else
//            {
//                dp[i].push_back(min(left, right));
//            }
//
//            if (i == n - 1 && j == 0) {
//                minnum = dp[i][j];
//            }
//            if (i == n - 1) {
//                minnum = min(minnum, dp[i][j]);
//            }
//        }
//    }
//
//    return minnum;
//}

//int numTrees(int n) {
//    if (n < 2)
//    {
//        return n;
//    }
//
//    vector<int> dp(n + 1, 0);
//    dp[0] = 1;
//    dp[1] = 1;
//
//    for (int i = 2; i <= n; i++)
//        for (int j = 1; j <= i; j++)
//            dp[i] += dp[j - 1] * dp[i - j];
//
//    return dp[n];
//}

//int searchInsert(vector<int>& nums, int target) {
//    int n = nums.size();
//    if (n == 0)
//    {
//        return 0;
//    }
//
//    int bIndex = 0;
//    int eIndex = n - 1;
//
//    /*while (bIndex < eIndex)
//    {
//        int mid = (bIndex + eIndex) / 2;
//        if (target > nums[mid])
//        {
//            bIndex = mid + 1;
//        }
//        else if (target == nums[mid])
//        {
//            return mid;
//        }
//        else 
//        {
//            eIndex = mid;
//        }
//    }
//
//    if (nums[bIndex] < target)
//    {
//        return n;
//    }
//
//    return bIndex;*/
//
//    while (bIndex <= eIndex)
//    {
//        if (target == nums[bIndex])
//        {
//            return bIndex;
//        }
//
//        if (target == nums[eIndex])
//        {
//            return eIndex;
//        }
//
//        if (target > nums[bIndex] && target < nums[eIndex])
//        {
//            ++bIndex;
//            --eIndex;
//            continue;
//        }
//
//        if (target < nums[bIndex])
//        {
//            int r = bIndex;
//            return r;
//        }
//
//        if (target > nums[eIndex])
//        {
//            int r = eIndex + 1;
//            return r;
//        }
//    }
//
//    return bIndex;
//}

//void rotate(vector<vector<int>>& matrix) {
//    int n = matrix.size();
//    int c = 0;
//    if (n % 2 == 0)
//    {
//        c = n - 1;
//    }
//    else
//    {
//        c = n - c / 2;
//    }
//
//    if (c < 0)
//    {
//        return;
//    }
//
//    for (int i = 0; i < c; ++i)
//    {
//        int step = n - 1 - 3 * i;
//        if (step < 0)
//        {
//            step = 0;
//        }
//
//        for (int j = 0; j < n; ++j)
//        {
//            for (int k = 0; k < n; ++k)
//            {
//                matrix[j][k]
//            }
//        }
//    }
//}



 
//bool isSameTree(TreeNode* p, TreeNode* q) {
//    if (p == nullptr)
//        return q == nullptr;
//    if (q == nullptr)
//        return p == nullptr;
//
//    if (p->val != q->val)
//    {
//        return false;
//    }
//
//    return isSameTree(p->left,q->left) && isSameTree(p->right, q->right);
//}

//bool isValidBST(TreeNode* root) {
//
//    stack<TreeNode*> stack;
//    long long inorder = (long long)INT_MIN - 1;
//
//    while (!stack.empty() || root != nullptr) {
//        while (root != nullptr) {
//            stack.push(root);
//            root = root->left;
//        }
//        root = stack.top();
//        stack.pop();
//
//        if (root->val <= inorder) return false;
//        inorder = root->val;
//        root = root->right;
//    }
//    return true;
//}

//bool check(TreeNode* left, TreeNode* right)
//{
//    if (left == nullptr)
//        return right == nullptr;
//    if (right == nullptr)
//        return left == nullptr;
//
//    if (left->val != right->val)
//    {
//        return false;
//    }
//
//    return left->val == right->val && check(left->left, right->right) && check(left->right, right->left);
//}
//
//bool isSymmetric(TreeNode* root) {
//
//    return check(root->left, root->right);
//}

//vector<vector<int>> levelOrder(TreeNode* root) {
//
//    vector<vector<int>> result;
//    if (root == nullptr)
//    {
//        return result;
//    }
//
//    queue<TreeNode*> roots;
//    roots.push(root);
//    while (!roots.empty())
//    {
//        int currentLevelSize = roots.size();
//        result.push_back(vector <int>());
//        for (int i = 0; i < currentLevelSize; ++i) {
//            auto node = roots.front(); roots.pop();
//            result.back().push_back(node->val);
//            if (node->left) roots.push(node->left);
//            if (node->right) roots.push(node->right);
//        }
//    }
//
//    return result;
//}

// 中序遍历循环写法
//vector<int> inorderTraversal(TreeNode* root) {
//    vector<int> result;
//    if (root == nullptr)
//    {
//        return result;
//    }
//
//    stack<TreeNode*> roots;
//
//    while (!roots.empty() || root != nullptr) {
//        while (root != nullptr)
//        {
//            roots.push(root);
//            root = root->left;
//        }
//
//        root = roots.top(); 
//        roots.pop();
//
//        if (root) {
//            result.push_back(root->val);
//            root = root->right;
//        }
//    }
//
//    return result;
//}

// 中序遍历递归写法
//vector<int> result;
//vector<int> inorderTraversal(TreeNode* root) {
//    
//    if (root != nullptr)
//    {
//        inorderTraversal(root->left);
//        result.push_back(root->val);
//        inorderTraversal(root->right);
//    }
//
//    return result;
//}


//int maxDepth(TreeNode* root) {
//    int depth = 0;
//    if (root != nullptr)
//    {
//        depth+=max(maxDepth(root->left), maxDepth(root->right));
//    }
//
//    return depth;
//}

//TreeNode* buildTree(vector<int>& preorder, vector<int>& inorder) {
//    int m = preorder.size();
//    if (m == 0) 
//        return NULL;
//    TreeNode* root = new TreeNode(preorder[0]);
//    vector<int> preorder_left, inorder_left, preorder_right, inorder_right;
//
//    int i;
//    //构造左子树的中序遍历
//    for (i = 0; i < m; i++) {
//        if (inorder[i] == root->val) break;
//        inorder_left.push_back(inorder[i]);
//    }
//    //构造右子树的中序遍历
//    for (i = i + 1; i < m; i++) {
//        inorder_right.push_back(inorder[i]);
//    }
//
//    int n = inorder_left.size();
//    for (int j = 1; j < m; j++) {
//        //构造左子树的前序遍历
//        if (j <= n)
//            preorder_left.push_back(preorder[j]);
//        //构造右子树的前序遍历
//        else preorder_right.push_back(preorder[j]);
//    }
//    root->left = buildTree(preorder_left, inorder_left);
//    root->right = buildTree(preorder_right, inorder_right);
//    return root;
//}

//int maxCoins(vector<int>& nums) {
//    int n = nums.size();
//    vector<int> tmp(n + 2);
//    tmp[0] = tmp[n + 1] = 1;
//    for (int i = 1; i <= n; i++) tmp[i] = nums[i - 1];
//    vector<vector<int>> dp(n + 2, vector<int>(n + 2, 0));
//    for (int step = 1; step <= n; step++) {
//        for (int i = 1; i + step - 1 <= n; i++) {
//            int j = i + step - 1;
//            for (int k = i; k <= j; k++) {
//                dp[i][j] = max(dp[i][j], dp[i][k - 1] + dp[k + 1][j] + tmp[i - 1] * tmp[k] * tmp[j + 1]);
//            }
//        }
//    }
//    return dp[1][n];
//}

//vector<int> twoSum(vector<int>& numbers, int target) {
//    int low = 0, high = numbers.size() - 1;
//    while (low < high) 
//    {
//        int sum = numbers[low] + numbers[high];
//        if (sum == target) 
//        {
//            return { low + 1, high + 1 };
//        }
//        else if (sum < target) 
//        {
//            ++low;
//        }
//        else 
//        {
//            --high;
//        }
//    }
//    return { -1, -1 };
//}

//vector<TreeNode*> buildTree(int begin, int end)
//{
//    if (begin > end)
//    {
//        return { nullptr };
//    }
//
//    vector<TreeNode*> allTrees;
//
//    for (int i = begin; i <= end; ++i)
//    {
//        vector<TreeNode*> leftTrees = buildTree(begin, i - 1);
//        vector<TreeNode*> rightTrees = buildTree(i+1, end);
//
//        for (auto& left : leftTrees)
//        {
//            for (auto& right : rightTrees)
//            {
//                TreeNode* root = new TreeNode(i);
//                root->left = left;
//                root->right = right;
//                allTrees.emplace_back(root);
//            }
//        }
//    }
//
//    return allTrees;
//}
//
//vector<TreeNode*> generateTrees(int n) {
//    if (n == 0)
//    {
//        return {};
//    }
//
//    return buildTree(1, n);
//}

//TreeNode* copyTree(TreeNode* tree)
//{
//    if (tree == nullptr) return nullptr;
//    TreeNode* newTree = new TreeNode(tree->val);
//    newTree->left = copyTree(tree->left);
//    newTree->right = copyTree(tree->right);
//
//    return newTree;
//}
//
//vector<TreeNode*> generateTrees(int n)
//{
//    if (n == 0) return {};
//    vector<TreeNode*> preTree;
//    preTree.emplace_back(nullptr);
//    for (int i = 1; i <= n; ++i)
//    {
//        vector<TreeNode*> currTree;
//        for (auto& root : preTree)
//        {
//            {
//                // 当根节点插入
//                TreeNode* newTree = new TreeNode(i);
//                newTree->left = root;
//                currTree.emplace_back(newTree);
//            }
//            
//            {
//                // 当右节点插入
//                for (int j = 0; j < n; ++j)
//                {
//                    TreeNode* tempTree = copyTree(root);
//                    TreeNode* right = tempTree;
//
//                    // 根据j，逐个寻找右节点
//                    for (int k = 0; k < j; ++k)
//                    {
//                        if (right == nullptr)
//                        {
//                            break;
//                        }
//
//                        right = right->right;
//                    }
//
//                    if (right == nullptr)
//                    {
//                        break;
//                    }
//
//                    // 替换原来的右节点，之前右节点都是新节点的左节点
//                    TreeNode* rightTree = right->right;
//                    TreeNode* newTree = new TreeNode(i);
//                    right->right = newTree;
//                    newTree->left = rightTree;
//                    currTree.emplace_back(tempTree);
//                }
//            }
//        }
//
//        preTree = currTree;
//    }
//
//    return preTree;
//}

//int minArray(vector<int>& numbers) {
//
//    int n = numbers.size();
//    if (n == 0)
//    {
//        return 0;
//    }
//    if (n == 1)
//    {
//        return numbers[0];
//    }
//
//    int b = 0;
//    int e = n - 1;
//    int mid = (e + b) / 2;
//    int minNum = numbers[0];
//
//    while (b != mid && e != mid)
//    {
//        if (numbers[b] > numbers[mid]) 
//        {
//            minNum = min(minNum, numbers[mid]);
//            e = mid;
//        }
//        else if(numbers[b] < numbers[mid])
//        {
//            minNum = min(minNum, numbers[b]);
//            b = mid;
//        }
//        else
//        {
//            ++b;
//        }
//
//        mid = (e + b) / 2;
//    }
//
//    minNum = min(minNum, min(numbers[b], numbers[e]));
//
//    return minNum;
//}

//int minArray(vector<int>& numbers) 
//{
//    int result = numbers[0];
//    int n = numbers.size();
//    for (int i = 1; i < n; ++i)
//    {
//        if (numbers[i] < result)
//        {
//            return numbers[i];
//        }
//    }
//
//    return result;
//}

//int minPathSum(vector<vector<int>>& grid) {
//    int n = grid.size();
//    int m = grid[0].size();
//
//    vector<vector<int>> dp(n, vector<int>(m,0));
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = 0; j < m; ++j)
//        {
//            dp[i][j] = grid[i][j];
//            if (i - 1 >= 0 && j - 1 >= 0)
//            {
//                dp[i][j] += min(dp[i - 1][j], dp[i][j - 1]);
//            }
//            else if (i - 1 >= 0)
//            {
//                dp[i][j] += dp[i - 1][j];
//            }
//            else if (j - 1 >= 0)
//            {
//                dp[i][j] += dp[i][j - 1];
//            }
//        }
//    }
//
//    return dp[n - 1][m - 1];
//}

//bool divisorGame(int N) {
//    return N % 2 == 0;
//}

//int helper(vector<int>& nums, int begin, int end, int d)
//{
//    if (begin > end)
//    {
//        return 0;
//    }
//
//    int result;
//    if (--d <= 0)
//    {
//        result = 0;
//        for (int i = begin; i < end; ++i)
//        {
//            result += nums[i];
//        }
//        return result;
//    }
//
//    result = INT_MAX;
//    for (int i = begin; end - (i + 1) >= d - 1; ++i)
//    {
//        int maxTemp = 0;
//        maxTemp = max(helper(nums, begin, i + 1, 0), helper(nums, i + 1, end, d));
//        result = min(result, maxTemp);
//    }
//
//    return result;
//}
//
//int splitArray(vector<int>& nums, int m) {
//    int n = nums.size();
//    int result = 0;
//    if (m == 1)
//    {
//        for (int i = 0; i < n; ++i)
//        {
//            result += nums[i];
//        }
//        return result;
//    }
//    if (n <= m)
//    {
//        for (int i = 0; i < n; ++i)
//        {
//            result = max(result, nums[i]);
//        }
//
//        return result;
//    }
//    result = INT_MAX;
//    for (int i = 0; n - (i + 1) >= m - 1; ++i)
//    {
//        int maxTemp = 0;
//        maxTemp = max(helper(nums, 0, i+1, 0), helper(nums, i + 1, n, m - 1));
//        result = min(result, maxTemp);
//    }
//
//    return result;
//}

//int splitArray(vector<int>& nums, int m) {
//    long l = nums[0], h = 0;//int类型在这里不合适，因为h可能会超过int类型能表示的最大值
//    for (auto i : nums)
//    {
//        h += i;
//        l = l > i ? l : i;
//    }
//    while (l < h)
//    {
//        long mid = (l + h) / 2;
//        long temp = 0;
//        int cnt = 1;//初始值必须为1
//        for (auto i : nums)
//        {
//            temp += i;
//            if (temp > mid)
//            {
//                temp = i;
//                ++cnt;
//            }
//        }
//        if (cnt > m)
//            l = mid + 1;
//        else
//            h = mid;
//    }
//    return l;
//}

// 只能买卖一次
//int maxProfit(vector<int>& prices) {
//    // 暴力解法
//    /*int n = prices.size();
//
//    int maxVaule = 0;
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = i + 1; j < n; ++j)
//        {
//            maxVaule = max(prices[j] - prices[i], maxVaule);
//        }
//    }
//
//    return maxVaule;*/
//
//    // 动态规划
//    //int n = prices.size();
//    //if (n == 0)
//    //{
//    //    return 0;
//    //}
//
//    //int minV = 0;
//    //int maxV = -prices[0];
//    //for (int i = 1; i < n; ++i)
//    //{
//    //    // 卖出利润
//    //    minV = max(minV, maxV + prices[i]);
//    //    // 买入利润
//    //    maxV = max(maxV, -prices[i]);
//    //}
//
//    //return max(minV,maxV);
//
//    // 贪心
//    int n = prices.size();
//    if (n == 0)
//    {
//        return 0;
//    }
//
//    int minV = INT_MAX;
//    int result = 0;
//    for (int i = 0; i < n; ++i)
//    {
//        if (prices[i] < minV)
//        {
//            minV = prices[i];
//        }
//        int temp = prices[i] - minV;
//        if (temp > result)
//        {
//            result = temp;
//        }
//    }
//
//    return result;
//}

// 可以多次买卖
//int maxProfit(vector<int>& prices) {
//	// 贪心
//	/*int n = prices.size();
//	if (n == 0)
//	{
//		return 0;
//	}
//
//    int result = 0;
//    for (int i = 1; i < n; ++i)
//    {
//        result += max(0, prices[i] - prices[i-1]);
//    }
//
//	return result;*/
//
//	// 动态规划
//	int n = prices.size();
//	if (n == 0)
//	{
//		return 0;
//	}
//
//	int sell = 0;
//	int buy = -prices[0];
//	for (int i = 1; i < n; ++i)
//	{
//		// 卖出利润
//        sell = max(sell, buy + prices[i]);
//		// 买入利润
//        buy = max(buy, sell -prices[i]);
//	}
//
//	return max(sell, buy);
//}

//bool canJump(vector<int>& nums) {
//    int n = nums.size();
//    if (n == 1)
//    {
//        return true;
//    }
//
//    if (nums[0] == 0)
//    {
//        return false;
//    }
//
//    int maxIndex = 0;
//    for (int i = 0; i < n; ++i)
//    {
//        if (i <= maxIndex)
//        {
//            maxIndex = max(maxIndex, i + nums[i]);
//        }
//        if (maxIndex >= n - 1)
//        {
//            return true;
//        }
//    }
//
//    return false;
//}

//int jump(vector<int>& nums) {
//	int n = nums.size();
//	if (n == 1)
//	{
//		return 0;
//	}
//
//    int count = 0;
//    int maxIndex = 0;
//
//    for (int i = 0; i < n; ++i)
//    {
//        int maxT = i + nums[i];
//        for (int j = i + 1; j <= maxIndex; ++j)
//        {
//            maxT = max(nums[j] + j, maxT);
//        }
//        maxIndex = maxT;
//        
//        ++count;
//        if (maxIndex >= n - 1)
//        {
//            break;
//        }
//    }
//
//    return count;
//}

//int jump(vector<int>& nums) {
//    int n = nums.size();
//    if (n == 1)
//    {
//        return 0;
//    }
//
//    int end = 0;
//    int count = 0;
//    int maxIndex = 0;
//
//    for (int i = 0; i < n - 1; ++i)
//    {
//        maxIndex = max(maxIndex, nums[i] + i);
//        if (i == end)
//        {
//            end = maxIndex;
//            ++count;
//        }
//    }
//
//    return count;
//}

//vector<int> divingBoard(int shorter, int longer, int k) {
//    if (k == 0)
//    {
//        return {};
//    }
//
//    if (shorter == longer)
//    {
//        return { shorter * k };
//    }
//    vector<int> result(k + 1, 0);
//    for (int i = 0; i <= k; ++i)
//    {
//        result[i] = (shorter * (k - i) + longer * i);
//    }
//
//    return result;
//}


//struct TrieTree {
//    bool isWord;
//    vector<TrieTree*> nextNode;
//
//    TrieTree() { nextNode.resize(26); isWord = false; }
//};
//
//void insertTrieNode(TrieTree* root, string str) {
//    if (root == nullptr || str.size() == 0)
//    {
//        return;
//    }
//    TrieTree* curr = root;
//    for (int i = str.length()-1; i >=0 ; --i) {
//        //如果该分支不存在，创建一个新节点
//        int c = str[i] - 'a';
//        if (curr->nextNode[c] == nullptr)
//        {
//            curr->nextNode[c] = new TrieTree();
//        }
//        curr = curr->nextNode[c];
//    }
//    curr->isWord = true;
//}
//
//
//vector<int> searchTrie(TrieTree* tree, string str,int endPos) {
//    vector<int> indices;
//    TrieTree* curr = tree;
//    for (int i = endPos; i>=0 ; --i)
//    {
//        int c = str[i] - 'a';
//        if (curr->nextNode[c] == nullptr)
//        {
//            break;
//        }
//        curr = curr->nextNode[c];
//        if (curr->isWord)
//        {
//            indices.emplace_back(i);
//        }
//    }
//    return indices;
//}
//
//
//int respace(vector<string>& dictionary, string sentence) {
//
//    TrieTree* root = new TrieTree();
//    int d_size = dictionary.size();
//    for (int i = 0; i < d_size; ++i)
//    {
//        insertTrieNode(root, dictionary[i]);
//    }
//
//    int sen_len = sentence.size();
//    vector<int> dp(sen_len + 1);      //dp[i]表示sentence前i个字符里面未识别的字符数
//
//    dp[0] = 0;
//    for (int i = 1; i <= sen_len; i++)
//    {
//        dp[i] = dp[i - 1] + 1;
//        for (auto t : searchTrie(root, sentence, i - 1))
//        {
//            dp[i] = min(dp[i], dp[t]);
//        }
//    }
//
//    return dp[sen_len];
//}

//bool isSubsequence(string s, string t) {
//    if (s.size() == 0)
//    {
//        return true;
//    }
//
//    int count = 0;
//    for (auto c : t)
//    {
//        if (s[count] == c)
//        {
//            ++count;
//        }
//
//        if (count == s.size())
//        {
//            return;
//        }
//    }
//
//    return false;
//}

//int canCompleteCircuit(vector<int>& gas, vector<int>& cost) {
//	int n = gas.size();
//	for (int i = 0; i < n; ++i)
//	{
//		if (gas[i] >= cost[i])
//		{
//			int cost_all = 0;
//			for (int j = 0; j < n; ++j)
//			{
//				int index = (j + i) >= n ? ((j + i) - n) : i + j;
//				cost_all += gas[index] - cost[index];
//				if (cost_all < 0)
//				{
//					break;
//				}
//			}
//
//			if (cost_all >= 0)
//			{
//				return i;
//			}
//		}
//	}
//
//	return -1;
//}

//int canCompleteCircuit(vector<int>& gas, vector<int>& cost) {
//    int curstart = 0;
//    int totalcost = 0;
//    int curcost = 0;
//    int size = gas.size();
//    for (int i = 0; i < size; i++)
//    {
//        totalcost = totalcost + gas[i] - cost[i];
//        curcost = curcost + gas[i] - cost[i];
//        if (curcost < 0)
//        {
//            curstart = i + 1;
//            curcost = 0;
//        }
//    }
//    return totalcost >= 0 ? curstart % size : -1;
//}

//int findContentChildren(vector<int>& g, vector<int>& s) {
//    int sn = s.size();
//    if(sn == 0)
//    {
//        return 0;
//    }
//    sort(g.begin(), g.end());
//    sort(s.begin(), s.end());
//    int j = sn - 1;
//    int count = 0;
//    for (int i = g.size() - 1; i >= 0; --i)
//    {
//        if (s[j] >= g[i]) 
//        {
//            --j;
//            ++count;
//        }
//
//        if (j < 0)
//        {
//            j = 0;
//            break;
//        }
//    }
//
//    return count;
//}

//int candy(vector<int>& ratings) {
//    int n = ratings.size();
//    if (n == 0)
//    {
//        return n;
//    }
//
//    vector<int> dp(n, 1);
//    for (int i = 1; i < n-1; ++i)
//    {
//        if (ratings[i] > ratings[i - 1] && ratings[i] > ratings[i + 1])
//        {
//            dp[i] += max(dp[i - 1], dp[i + 1]);
//        }
//        else if (ratings[i] > ratings[i - 1])
//        {
//            dp[i] += dp[i - 1];
//        }
//        else if (ratings[i] > ratings[i + 1])
//        {
//            dp[i] += dp[i + 1];
//        }
//        else if (ratings[i] < ratings[i - 1] && ratings[i] < ratings[i + 1])
//        {
//            ++dp[i - 1];
//            ++dp[i + 1];
//        }
//    }
//
//    int result = 0;
//    for(auto i : dp)
//    {
//        result += i;
//    }
//
//
//    return result;
//}

//vector<string> generateParenthesis(int n) {
//    vector<string> res;
//    int lc = 0, rc = 0;
//    dfs(res, "", n, lc, rc);
//    return res;
//}
//void dfs(vector<string>& res, string path, int n, int lc, int rc) {
//    if (rc > lc || lc > n || rc > n) return;
//    if (lc == rc && lc == n) {
//        res.push_back(path);
//        return;
//    }
//    dfs(res, path + '(', n, lc + 1, rc);
//    dfs(res, path + ')', n, lc, rc + 1);
//}

//void setBlock(vector<vector<int>>& temp,int x,int y, int block)
//{
//    int n = temp.size();
//    for (int j = x + 1; j < n; ++j)
//    {
//        for (int k = 0; k < n; ++k)
//        {
//            if (temp[j][k] == 1)
//            {
//                continue;
//            }
//            if (temp[j][k] == 0 && (k == y || (abs(j-x) == abs(k - y))))
//            {
//                temp[j][k] = block;
//            }
//        }
//    }
//}
//
//void cleanBlock(vector<vector<int>>& temp, int upper, int block)
//{
//    int n = temp.size();
//    for (int j = upper; j < n; ++j)
//    {
//        for (int k = 0; k < n; ++k)
//        {
//            if (block == temp[j][k])
//            {
//                temp[j][k] = 0;
//            }
//        }
//    }
//}
//
//void buildString(vector<vector<string>>& result,vector<vector<int>> &temp, int upper,int n)
//{
//    if (upper == n)
//    {
//        vector<string> r;
//        for (int j = 0; j < n; ++j)
//        {
//            int block_count = 0;
//            string t(n,'.');
//            for (int k = 0; k < n; ++k)
//            {
//                if (temp[j][k] >= 2)
//                {
//                    ++block_count;
//                }
//                else if (temp[j][k] == 1)
//                {
//                    t[k] = 'Q';
//                }
//            }
//            if (block_count == n)
//            {
//                return;
//            }
//
//            r.push_back(t);
//        }
//
//        result.push_back(r);
//        return;
//    }
//
//    for (int i = 0; i < n; ++i)
//    {
//        if (temp[upper][i] != 0)
//        {
//            continue;
//        }
//
//        temp[upper][i] = 1;
//        setBlock(temp, upper, i, upper + 2);
//        buildString(result, temp, upper + 1, n);
//        temp[upper][i] = 0;
//        cleanBlock(temp, upper + 1, upper + 2);
//    }
//}
//
//vector<vector<string>> solveNQueens(int n) {
//
//    if (n == 0)
//    {
//        return { {} };
//    }
//
//    vector<vector<string>> result;
//    for (int i = 0; i < n; ++i)
//    {
//        vector<vector<int>> temp(n, vector<int>(n, 0));
//        temp[0][i] = 1;
//        setBlock(temp, 0, i, 2);
//        buildString(result, temp,1,n);
//    }
//
//    return result;
//}

//int distributeCandies(vector<int>& candies) {
//    set<int> m;
//    for (auto c : candies)
//    {
//        m.insert(c);
//    }
//    
//    return min(m.size(), (candies.size() / 2));
//}

//vector<int> intersect(vector<int>& nums1, vector<int>& nums2) {
//    
//    map<int, int> m1;
//    for (auto n : nums1)
//    {
//        if (m1.find(n) != m1.end())
//        {
//            ++m1[n];
//        }
//        else
//        {
//            m1[n] = 1;
//        }
//    }
//
//    vector<int> result;
//    for (auto n : nums2)
//    {
//        auto f = m1.find(n);
//        if (f != m1.end())
//        {
//            result.push_back(n);
//            --f->second;
//            if (f->second == 0)
//            {
//                m1.erase(f);
//            }
//        }
//    }
//
//    return result;
//}

//vector<vector<int>> pairSums(vector<int>& nums, int target) {
//    unordered_map<int, int> m;
//    vector<vector<int>> result;
//    for (auto n : nums)
//    {
//        if (m.find(n) != m.end())
//        {
//            ++m[n];
//        }
//        else
//        {
//            m[n] = 1;
//        }
//    }
//
//    for (auto t : m)
//    {
//        if (t.second == 0)
//        {
//            continue;
//        }
//        int d = target - t.first;
//        auto tt = m.find(d);
//        int times = 0;
//        vector<int> r;
//        if (tt != m.end() && tt->second != 0)
//        {
//            if (tt->first == t.first)
//            {
//                times = t.second / 2;
//            }
//            else
//            {
//                times = min(t.second, tt->second);
//                tt->second -= times;
//            }
//            t.second -= times;
//            r.push_back(t.first);
//            r.push_back(tt->first);
//        }
//
//        for (int i = 0; i < times; ++i)
//        {
//            result.push_back(r);
//        }
//    }
//
//    return result;
//}

//int islandPerimeter(vector<vector<int>>& grid) {
//    int n = grid.size();
//    unordered_map<int, int> m;
//    for (int i = 0; i < n; ++i)
//    {
//        int count = 0;
//        int nn = grid[i].size();
//        for (int j = 0; j < nn; ++j)
//        {
//            if (grid[i][j] == 1)
//            {
//                ++count;
//            }
//        }
//
//        m[i] = count;
//    }
//
//    int result = (m[0] > 0) ? m[0] * 4 - ((m[0] - 1) * 2) : 0;
//    for (int i = 1; i < n; ++i)
//    {
//        if (m[i] == 0)
//        {
//            break;
//        }
//        int dis = min(m[i - 1], m[i]);
//        result += m[i] * 4 - ((m[i] - 1) * 2) - dis * 2;
//    }
//
//    return result;
//}

//void reverseString(vector<char>& s) {
//    int n = s.size();
//    if (n <= 1)
//    {
//        return;
//    }
//
//    int l = 0;
//    int r = n - 1;
//    while (l < r)
//    {
//        char temp = s[l];
//        s[l] = s[r];
//        s[r] = temp;
//
//        ++l;
//        --r;
//    }
//}

//vector<vector<int>> threeSum(vector<int>& nums) {
//    int n = nums.size();
//    sort(nums.begin(), nums.end());
//    vector<vector<int>> ans;
//    // 枚举 a
//    for (int first = 0; first < n; ++first) {
//        // 需要和上一次枚举的数不相同
//        if (first > 0 && nums[first] == nums[first - 1]) {
//            continue;
//        }
//        // c 对应的指针初始指向数组的最右端
//        int third = n - 1;
//        int target = -nums[first];
//        // 枚举 b
//        for (int second = first + 1; second < n; ++second) {
//            // 需要和上一次枚举的数不相同
//            if (second > first + 1 && nums[second] == nums[second - 1]) {
//                continue;
//            }
//            // 需要保证 b 的指针在 c 的指针的左侧
//            while (second < third && nums[second] + nums[third] > target) {
//                --third;
//            }
//            // 如果指针重合，随着 b 后续的增加
//            // 就不会有满足 a+b+c=0 并且 b<c 的 c 了，可以退出循环
//            if (second == third) {
//                break;
//            }
//            if (nums[second] + nums[third] == target) {
//                ans.push_back({ nums[first], nums[second], nums[third] });
//            }
//        }
//    }
//    return ans;
//}

//string addStrings(string num1, string num2) {
//    string result;
//
//    int l1 = num1.size() -1;
//    int l2 = num2.size() -1;
//    bool up = false;
//    while (l1 >= 0 || l2 >= 0)
//    {
//        int t1 = 0;
//        if (l1 >= 0)
//        {
//            t1 = num1[l1] - 48;
//        }
//        
//        int t2 = 0;
//        if (l2 >= 0)
//        {
//            t2 = num2[l2] - 48;
//        }
//        int tr = t1 + t2 + ((up) ? 1 : 0);
//        if (tr >= 10)
//        {
//            up = true;
//            tr = tr % 10;
//        }
//        else
//        {
//            up = false;
//        }
//
//        char tc = tr + 48;
//        result.push_back(tc);
//        --l1;
//        --l2;
//    }
//
//    if (up)
//    {
//        result.push_back('1');
//    }
//
//    reverse(result.begin(), result.end());
//
//    return result;
//}

//vector<vector<int>> edges;
//vector<int> visited;
//bool valid = true;
//
//void dfs(int u)
//{
//    visited[u] = 1;
//    for (auto& i : edges[u])
//    {
//        // 没走过就向下找
//        if (visited[i] == 0)
//        {
//            dfs(i);
//            if (!valid)
//            {
//                return;
//            }
//        }
//        // 已经走过，直接返回错误
//        else if(visited[i] == 1)
//        {
//            valid = false;
//            return;
//        }
//    }
//
//    visited[u] = 2;
//}
//
//bool canFinish(int numCourses, vector<vector<int>>& prerequisites) {
//    
//    edges.resize(numCourses);
//    visited.resize(numCourses);
//
//    for (auto& info : prerequisites)
//    {
//        edges[info[1]].push_back(info[0]);
//    }
//
//    for (int i = 0; i < numCourses && valid; ++i)
//    {
//        if (!visited[i])
//        {
//            dfs(i);
//        }
//    }
//
//    return valid;
//}

//vector<vector<int>> edges;
//vector<int> visited;
//vector<int> result;
//bool valid = true;
//
//void dfs(int u)
//{
//    visited[u] = 1;
//    result.push_back(u);
//    for (auto& i : edges[u])
//    {
//        
//        // 没走过就向下找
//        if (visited[i] == 0)
//        {
//            dfs(i);
//            if (!valid)
//            {
//                return;
//            }
//        }
//        // 已经走过，直接返回错误
//        else if (visited[i] == 1)
//        {
//            valid = false;
//            return;
//        }
//    }
//    
//    visited[u] = 2;
//}
//
//vector<int> findOrder(int numCourses, vector<vector<int>>& prerequisites) {
//
//    edges.resize(numCourses);
//    visited.resize(numCourses);
//
//    for (auto& info : prerequisites)
//    {
//        edges[info[1]].push_back(info[0]);
//    }
//
//    for (int i = 0; i < numCourses && valid; ++i)
//    {
//        if (!visited[i])
//        {
//            dfs(i);
//        }
//    }
//
//    if (!valid)
//    {
//        return { {} };
//    }
//
//    reverse(result.begin(), result.end());
//    return result;
//}

//unordered_map<TreeNode*, int> d,n;
//
//void dfs(TreeNode* root)
//{
//    if (root == nullptr)
//    {
//        return;
//    }
//
//    dfs(root->left);
//    dfs(root->right);
//
//    // 偷这个：当前值加上不偷左右两个节点的值
//    d[root] = root->val + n[root->left] + n[root->right];
//
//    // 不偷这个：取偷或不偷左右节点时的最大值
//    n[root] = max(d[root->left], n[root->left])+ max(d[root->right], n[root->right]);
//}
//
//int rob(TreeNode* root) {
//    if (root == nullptr)
//    {
//        return 0;
//    }
//    dfs(root);
//    
//    int result = max(d[root], n[root]);
//    return result;
//}

//unordered_map<string, int> strings;
//
//int findWord(const string& s, int left, int right) {
//    auto iter = strings.find(s.substr(left, right - left + 1));
//    return iter == strings.end() ? -1 : iter->second;
//}
//
//bool isPalindrome(const string& s, int left, int right) {
//    int len = right - left + 1;
//    for (int i = 0; i < len / 2; i++) {
//        if (s[left + i] != s[right - i]) {
//            return false;
//        }
//    }
//    return true;
//}
//
//vector<vector<int>> palindromePairs(vector<string>& words) {
//    int n = words.size();
//
//    // 存倒串
//    for (int i = 0; i < n; ++i)
//    {
//        string temp = words[i];
//        reverse(temp.begin(), temp.end());
//        strings[temp] = i;
//    }
//
//    vector<vector<int>> result;
//    for (int i = 0; i < n; ++i)
//    {
//        int m = words[i].size();
//        if (!m) {
//            continue;
//        }
//
//        for (int j = 0; j <= m; ++j)
//        {
//            //左边匹配列表中的其他串，即右边为回文串
//            if (isPalindrome(words[i], j, m - 1))
//            {
//                int index = findWord(words[i], 0, j - 1);
//                if (index != -1 && index != i)
//                {
//                    result.push_back({ i,index });
//                }
//            }
//
//            //右边匹配列表中的其他串，即左边为回文串
//            if (j && isPalindrome(words[i], 0, j - 1))
//            {
//                int index = findWord(words[i], j, m - 1);
//                if (index != -1 && index != i)
//                {
//                    result.push_back({ index,i });
//                }
//            }
//        }
//    }
//
//    return result;
//}

//void recoverTree(TreeNode* root) {
//    stack<TreeNode*> stack;
//
//    TreeNode* t1 = nullptr;
//    TreeNode* t2 = nullptr;
//    TreeNode* lower = nullptr;
//
//    while (!stack.empty() || root != nullptr)
//    {
//        while (root != nullptr)
//        {
//            stack.push(root);
//            root = root->left;
//        }
//
//        root = stack.top();
//        stack.pop();
//
//        if (lower != nullptr && root->val <= lower->val)
//        {
//            t2 = root;
//            if (t1 == nullptr)
//            {
//                t1 = lower;
//            }
//            else
//            {
//                break;
//            }
//        }
//
//        lower = root;
//        root = root->right;
//    }
//
//    swap(t1->val, t2->val);
//}

//bool isStart = false;
//TreeNode* startPos;
//bool subPathHelper(ListNode* origin, ListNode* head, TreeNode* root)
//{
//    if (head == nullptr)
//    {
//        return true;
//    }
//
//    if (root == nullptr)
//    {
//        return false;
//    }
//
//    ListNode* temp = head;
//    if (isStart)
//    {
//        if (root->val != head->val)
//        {
//            isStart = false;
//            return subPathHelper(origin, origin, startPos->left) || subPathHelper(origin, origin, startPos->right);
//        }
//        else
//        {
//            temp = head->next;
//        }
//    }
//    else
//    {
//        if (root->val == head->val)
//        {
//            temp = head->next;
//            startPos = root;
//            isStart = true;
//        }
//    }
//    return subPathHelper(origin, temp, root->left) || subPathHelper(origin, temp, root->right);
//}
//
//bool isSubPath(ListNode* head, TreeNode* root) {
//
//    isStart = false;
//    return subPathHelper(head, head, root);
//}

//bool flipEquiv(TreeNode* root1, TreeNode* root2) {
//    
//    if (root1 == root2)
//    {
//        return true;
//    }
//
//    if (root1 == nullptr || root2 == nullptr || root1->val != root2->val)
//    {
//        return false;
//    }
//
//    return flipEquiv(root1->left,root2->left) && flipEquiv(root1->right, root2->right)
//        || flipEquiv(root1->right, root2->left) && flipEquiv(root1->left, root2->right);
//}

//vector<int> postorder(Node* root) {
//
//    if (root == nullptr)
//    {
//        return {};
//    }
//
//    stack<Node*> stack;
//    vector<int> result;
//    stack.push(root);
//    while (!stack.empty())
//    {
//        Node* temp = stack.top();
//        stack.pop();
//        result.push_back(temp->val);
//        for (int i = 0; i < temp->children.size(); ++i)
//        {
//            stack.push(temp->children[i]);
//        }
//    }
//
//    reverse(result.begin(), result.end());
//    return result;
//}


//vector<string> result;
//vector<int> tip;
//
//void restoreIpHelper(string s, int segId, int segStart)
//{
//    if (segId == 4)
//    {
//        if (s.size() == segStart)
//        {
//            string ipAddr;
//            for (int i = 0; i < 4; ++i) {
//                ipAddr += to_string(tip[i]);
//                ipAddr += ".";
//            }
//            
//            ipAddr.erase(ipAddr.end() - 1);
//            result.push_back(move(ipAddr));
//        }
//        return;
//    }
//
//    // 字符串长度小于可分割数
//    if (s.size() == segStart)
//    {
//        return;
//    }
//
//    // 不能以0为开头，所以就单独为一位
//    if (s[segStart] == '0')
//    {
//        tip[segId] = 0;
//        restoreIpHelper(s, segId + 1, segStart + 1);
//    }
//
//    int addr = 0;
//    for (int i = segStart; i < s.size(); ++i)
//    {
//        addr = addr * 10 + (s[i] - '0');
//        if (addr > 0 && addr <= 0xff)
//        {
//            tip[segId] = addr;
//            restoreIpHelper(s, segId + 1, i + 1);
//        }
//        else
//        {
//            break;
//        }
//    }
//}
//
//vector<string> restoreIpAddresses(string s) {
//    if (s.empty() || s.size() < 4)
//    {
//        return {};
//    }
//    tip.resize(4);
//    restoreIpHelper(s, 0, 0);
//
//    return result;
//}


//int countBinarySubstrings(string s) {
//
//    int n = s.size();
//    int last = 0;
//    int curr = 0;
//    int result = 0;
//    char* c = &s[0];
//    for (int i = 0; i < n; ++i)
//    {
//        if (*c == s[i])
//        {
//            ++curr;
//        }
//        else
//        {
//            last = curr;
//            c = &s[i];
//            curr = 1;
//        }
//
//        if (last >= curr)
//        {
//            ++result;
//        }
//    }
//
//    return result;
//}

//int n, m;
//
//void dfs(vector<vector<char>>& board, int x, int y) {
//    if (x < 0 || x >= n || y < 0 || y >= m || board[x][y] != 'O') {
//        return;
//    }
//    board[x][y] = 'A';
//    if (x + 1 < n) dfs(board, x + 1, y);
//    if (x - 1 >= 0)dfs(board, x - 1, y);
//    if (y + 1 < m)dfs(board, x, y + 1);
//    if (y - 1 >= 0)dfs(board, x, y - 1);
//}
//
//void solve(vector<vector<char>>& board) {
//    n = board.size();
//    if (n == 0) {
//        return;
//    }
//    m = board[0].size();
//    for (int i = 0; i < n; i++) {
//        dfs(board, i, 0);
//        dfs(board, i, m - 1);
//    }
//    for (int i = 1; i < m - 1; i++) {
//        dfs(board, 0, i);
//        dfs(board, n - 1, i);
//    }
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < m; j++) {
//            if (board[i][j] == 'A') {
//                board[i][j] = 'O';
//            }
//            else if (board[i][j] == 'O') {
//                board[i][j] = 'X';
//            }
//        }
//    }
//}

//Node* cloneGraph(Node* node) {
//
//    if (node == nullptr)
//    {
//        return {};
//    }
//
//    Node* result = new Node(1);
//
//    if (node->neighbors.size() == 0)
//    {
//        return result;
//    }
//    
//    unordered_map<int, Node*> visited;
//    unordered_map<int, Node*> has;
//    queue<Node*> que;
//    que.push(node);
//    queue<Node*> temp;
//    temp.push(result);
//    has[result->val] = result;
//    while (!que.empty())
//    {
//        Node* root = que.front();
//        que.pop();
//
//        Node* n = temp.front();
//        temp.pop();
//
//        visited[root->val] = n;
//        for (auto& index : root->neighbors)
//        {
//            if (visited.find(index->val) == visited.end())
//            {
//                Node* t = has.find(index->val) != has.end() ? has.find(index->val)->second : new Node(index->val);
//                t->neighbors.push_back(n);
//                n->neighbors.push_back(t);
//                has[t->val] = t;
//
//                que.push(index);
//                temp.push(t);
//            }
//        }
//    }
//
//
//    return result;
//}

//string multiply(string num1, string num2)
//{
//    string result;
//    vector<int> str;
//    bool up = false;
//    int up_num = 0;
//    for (int i = num1.size() - 1 ; i >= 0; --i)
//    {
//        up = false;
//        up_num = 0;
//        str.clear();
//
//        int t1 = 0;
//        if (i >= 0)
//        {
//            t1 = num1[i] - 48;
//        }
//        for (int j = num2.size() - 1; j >= 0; --j)
//        {
//            int t2 = 0;
//            if (j >= 0)
//            {
//                t2 = num2[j] - 48;
//            }
//            int tr = t1 * t2 + up_num;
//            if (tr >= 10)
//            {
//                up = true;
//                up_num = tr / 10;
//                tr = tr % 10;
//            }
//            else
//            {
//                up = false;
//                up_num = 0;
//            }
//            str.push_back(tr);
//        }
//
//        if (up)
//        {
//            str.push_back(up_num);
//        }
//
//        reverse(str.begin(), str.end());
//
//        for (int k = 0; k < num1.size() - i - 1; ++k)
//        {
//            str.push_back(0);
//        }
//
//        string r_str;
//        int n1 = result.size() - 1;
//        int n2 = str.size() - 1;
//        up = false;
//        while (n1 >= 0 || n2 >= 0)
//        {
//            int t1 = 0;
//            if (n1 >= 0)
//            {
//                t1 = result[n1] - 48;
//            }
//
//            int t2 = 0;
//            if (n2 >= 0)
//            {
//                t2 = str[n2];
//            }
//            int tr = t1 + t2 + ((up) ? 1 : 0);
//            if (tr >= 10)
//            {
//                up = true;
//                tr = tr % 10;
//            }
//            else
//            {
//                up = false;
//            }
//
//            char tc = tr + 48;
//            r_str.push_back(tc);
//            --n1;
//            --n2;
//        }
//
//        if (up)
//        {
//            r_str.push_back('1');
//        }
//
//        reverse(r_str.begin(), r_str.end());
//        result = r_str;
//    }
//
//    if (result[0] == '0')
//    {
//        return "0";
//    }
//    return result;
//}

//const int dir[4][2] = { {0,1},{0,-1},{1,0},{-1,0} };
//void dfs(vector<vector<int>>& image, int x, int y, int color, int newColor) {
//    if (image[x][y] == color) {
//        image[x][y] = newColor;
//        for (int i = 0; i < 4; i++) {
//            int mx = x + dir[i][0], my = y + dir[i][1];
//            if (mx >= 0 && mx < image.size() && my >= 0 && my < image[0].size()) {
//                dfs(image, mx, my, color, newColor);
//            }
//        }
//    }
//}
//
//vector<vector<int>> floodFill(vector<vector<int>>& image, int sr, int sc, int newColor) {
//    int currColor = image[sr][sc];
//    if (currColor != newColor) {
//        dfs(image, sr, sc, currColor, newColor);
//    }
//    return image;
//}

//vector<vector<int>> floodFill(vector<vector<int>>& image, int sr, int sc, int newColor) {
//
//    int n = image.size();
//    int m = image[0].size();
//    int currColor = image[sr][sc];
//    if (currColor == newColor)
//    {
//        return image;
//    }
//
//    int dir[4][2] = { {0,1},{0,-1},{1,0},{-1,0} };
//    queue<pair<int, int>> que;
//    que.push({ sr,sc });
//    image[sr][sc] = newColor;
//    while (!que.empty())
//    {
//        int x = que.front().first;
//        int y = que.front().second;
//        que.pop();
//
//        for (int i = 0; i < 4; ++i)
//        {
//            int t_x = x + dir[i][0];
//            int t_y = y + dir[i][1];
//            if (t_x < 0 || t_y < 0 || t_x >= n || t_y >= m || image[t_x][t_y] == newColor || image[t_x][t_y] == 0)
//            {
//                continue;
//            }
//            image[t_x][t_y] = newColor;
//            que.push({ t_x,t_y });
//        }
//    }
//
//    return image;
//}

//TreeNode* sortedListToBST(ListNode* head) {
//
//}
//


//TreeNode* sortedArrayToBSTHelper(vector<int>& nums, int b, int e)
//{
//    if (b > e)
//    {
//        return nullptr;
//    }
//
//    int mid = (e + b) / 2;
//    TreeNode* root = new TreeNode(nums[mid]);
//    root->left = sortedArrayToBSTHelper(nums, b, mid - 1);
//    root->right = sortedArrayToBSTHelper(nums, mid + 1, e);
//    return root;
//}
//
//TreeNode* sortedArrayToBST(vector<int>& nums)
//{
//    return sortedArrayToBSTHelper(nums, 0, nums.size() - 1);
//}

//int countSubstrings(string s) {
//    int n = s.size();
//    int result = 0;
//    vector<vector<bool>> dp(n, vector<bool>(n, false));
//
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = 0; j <= i; ++j)
//        {
//            if (i == j)
//            {
//                dp[j][i] = true;
//                ++result;
//            }
//            else if (i - j == 1 && s[i] == s[j])
//            {
//                dp[j][i] = true;
//                ++result;
//            }
//            else if (i - j > 1 && s[i] == s[j] && dp[j + 1][i - 1])
//            {
//                dp[j][i] = true;
//                ++result;
//            }
//        }
//    }
//
//    return result;
//}

//int dir_x[8] = { 0, 1, 0, -1, 1, 1, -1, -1 };
//int dir_y[8] = { 1, 0, -1, 0, 1, -1, 1, -1 };
//
//void dfs(vector<vector<char>>& board, int x, int y) {
//    int cnt = 0;
//    for (int i = 0; i < 8; ++i) {
//        int tx = x + dir_x[i];
//        int ty = y + dir_y[i];
//        if (tx < 0 || tx >= board.size() || ty < 0 || ty >= board[0].size()) {
//            continue;
//        }
//        // 不用判断 M，因为如果有 M 的话游戏已经结束了
//        cnt += board[tx][ty] == 'M';
//    }
//    if (cnt > 0) {
//        // 规则 3
//        board[x][y] = cnt + '0';
//    }
//    else {
//        // 规则 2
//        board[x][y] = 'B';
//        for (int i = 0; i < 8; ++i) {
//            int tx = x + dir_x[i];
//            int ty = y + dir_y[i];
//            // 这里不需要在存在 B 的时候继续扩展，因为 B 之前被点击的时候已经被扩展过了
//            if (tx < 0 || tx >= board.size() || ty < 0 || ty >= board[0].size() || board[tx][ty] != 'E') {
//                continue;
//            }
//            dfs(board, tx, ty);
//        }
//    }
//}
//
//vector<vector<char>> updateBoard(vector<vector<char>>& board, vector<int>& click) {
//    int x = click[0], y = click[1];
//    if (board[x][y] == 'M') {
//        // 规则 1
//        board[x][y] = 'X';
//    }
//    else {
//        dfs(board, x, y);
//    }
//    return board;
//}

//int rangeBitwiseAnd(int m, int n) {
//    int result = 0;
//
//    if (m == n)
//    {
//        return m;
//    }
//
//    for (int i = 0; i < 31; ++i)
//    {
//        u_int temp1 = 1 << i;
//        u_int temp2 = 1 << (i + 1);
//        
//        if (m >= temp1 && n < temp2)
//        {
//            result = temp1 + rangeBitwiseAnd(m - temp1, n - temp1);
//            break;
//        }
//    }
//
//    return result;
//}

//vector<int> temp;
//vector<vector<int>> ans;
//
//void findSubsequencesDfs(int cur, int last, vector<int>& nums) {
//    if (cur == nums.size()) {
//        if (temp.size() >= 2) {
//            ans.push_back(temp);
//        }
//        return;
//    }
//    if (nums[cur] >= last) {
//        temp.push_back(nums[cur]);
//        findSubsequencesDfs(cur + 1, nums[cur], nums);
//        temp.pop_back();
//    }
//    if (nums[cur] != last) {
//        findSubsequencesDfs(cur + 1, last, nums);
//    }
//}
//
//vector<vector<int>> findSubsequences(vector<int>& nums) {
//    findSubsequencesDfs(0, INT_MIN, nums);
//    return ans;
//}

//vector<int> vis;
//int num;
//
//void dfs(vector<vector<int>>& rooms, int x) {
//    vis[x] = true;
//    num++;
//    for (auto& it : rooms[x]) {
//        if (!vis[it]) {
//            dfs(rooms, it);
//        }
//    }
//}
//
//bool canVisitAllRooms(vector<vector<int>>& rooms) {
//    int n = rooms.size();
//    num = 0;
//    vis.resize(n);
//    dfs(rooms, 0);
//    return num == n;
//}

//vector<vector<int>> levelOrderBottom(TreeNode* root) {
//    vector<vector<int>> result;
//    if (!root) {
//        return result;
//    }
//    queue<TreeNode*> q;
//    q.push(root);
//    vector<int> level;
//    while (!q.empty()) {
//
//        int size = q.size();
//        for (int i = 0; i < size; ++i) {
//            auto node = q.front();
//            q.pop();
//            level.push_back(node->val);
//            if (node->left) {
//                q.push(node->left);
//            }
//            if (node->right) {
//                q.push(node->right);
//            }
//        }
//        if (!level.empty()) result.push_back(level);
//        level.clear();
//    }
//    reverse(result.begin(), result.end());
//    return result;
//}

//class mycomparison {
//public:
//    bool operator()(const pair<int, int>& lhs, const pair<int, int>& rhs) {
//        return lhs.second > rhs.second;
//    }
//};
//
//    vector<int> topKFrequent(vector<int>& nums, int k) {
//        // 要统计元素出现频率
//        unordered_map<int, int> map; // map<nums[i],对应出现的次数>
//        for (int i = 0; i < nums.size(); i++) {
//            map[nums[i]]++;
//        }
//
//        // 对频率排序
//        // 定义一个小顶堆，大小为k
//        priority_queue<pair<int, int>, vector<pair<int, int>>, mycomparison> pri_que;
//        for (unordered_map<int, int>::iterator it = map.begin(); it != map.end(); it++) {
//            pri_que.push(*it);
//            if (pri_que.size() > k) {
//                pri_que.pop();
//            }
//        }
//
//        // 找出前K个高频元素，因为小顶堆先弹出的是最小的，所以倒叙来数值数组
//        vector<int> result(k);
//        for (int i = k - 1; i >= 0; i--) {
//            result[i] = pri_que.top().first;
//            pri_que.pop();
//        }
//        return result;
//
//    }

//vector<vector<int>> result;
//vector<int> temp;
//
//void combineHelper(int b, int e, int k)
//{
//    if (k == 0)
//    {
//        result.push_back(temp);
//        return;
//    }
//
//    for (int i = b; i <= e; ++i)
//    {
//        temp.push_back(i);
//        combineHelper(i + 1, e, k - 1);
//        temp.pop_back();
//    }
//}
//
//vector<vector<int>> combine(int n, int k) {
//
//    combineHelper(1, n, k);
//    return result;
//}

//vector<vector<int>> result;
//vector<int> temp;
//
//void combinationSumHelper(vector<int>& candidates, int b, int t)
//{
//    if (t == 0)
//    {
//        result.push_back(temp);
//        return;
//    }
//
//    int n = candidates.size();
//    for (int i = b; i < n; ++i)
//    {
//        int sub = t - candidates[i];
//        if (sub < 0)
//        {
//            continue;
//        }
//        temp.push_back(candidates[i]);
//        combinationSumHelper(candidates, i, sub);
//        temp.pop_back();
//    }
//}
//
//vector<vector<int>> combinationSum(vector<int>& candidates, int target) {
//
//
//    combinationSumHelper(candidates, 0, target);
//    return result;
//}

//vector<pair<int, int>> dir = { {0,1},{1,0}, {0,-1},{-1,0} };
//
//bool existHelper(vector<vector<char>>& board, string& word, vector<vector<bool>>& visited, int i,int j,int index) {
//    if (word[index] != board[i][j])
//    {
//        return false;
//    }
//    
//    if (index == word.size() - 1)
//    {
//        return true;
//    }
//
//    visited[i][j] = true;
//    for (int d = 0; d < 4; ++d)
//    {
//        int x = i + dir[d].first;
//        int y = i + dir[d].second;
//        if (x >= 0 && x < visited.size() && y >= 0 && y < visited[0].size() && !visited[x][y])
//        {
//            if (existHelper(board, word, visited, i, j, index + 1))
//            {
//                return true;
//            }
//        }
//    }
//
//    visited[i][j] = false;
//
//    return false;
//}
//
//bool exist(vector<vector<char>>& board, string word) {
//    
//    int n = board.size();
//    int m = board[0].size();
//    vector<vector<bool>> visited(n,vector<bool>(m,false));
//    for (int i = 0; i < n; ++i)
//    {
//        for (int j = 0; j < m; ++j)
//        {
//            if (existHelper(board, word, visited,i,j,0))
//            {
//                return true;
//            }
//        }
//    }
//    return false;
//}

//int removeDuplicates(vector<int>& nums) {
//    if (nums.size() == 0)
//    {
//        return 0;
//    }
//    int i = 0;
//    int j = 1;
//    int result = 1;
//    while (j < nums.size())
//    {
//        if (nums[i] != nums[j])
//        {
//            ++i;
//            swap(nums[i], nums[j]);
//            
//            ++result;
//        }
//        ++j;
//    }
//
//    return result;
//}

//int removeElement(vector<int>& nums, int val) {
//    int i = 0;
//    int n = nums.size();
//    while (i < n) {
//        if (nums[i] == val) {
//            nums[i] = nums[n - 1];
//            // reduce array size by one
//            n--;
//        }
//        else {
//            i++;
//        }
//    }
//    return n;
//}


//vector<int> preorderTraversal(TreeNode* root) {
//
//    if (root == nullptr)
//    {
//        return {};
//    }
//    vector<int> result;
//    stack<TreeNode*> stk;
//    stk.push(root);
//    while (!stk.empty() && root != nullptr)
//    {
//        TreeNode* root = stk.top();
//        stk.pop();
//
//        result.push_back(root->val);
//        if (root->right)
//        {
//            stk.push(root->right);
//        }
//        if (root->left)
//        {
//            stk.push(root->left);
//        }
//    }
//
//    return result;
//}

//vector<int> postorderTraversal(TreeNode* root) {
//    stack<TreeNode*> st;
//    vector<int> result;
//    st.push(root);
//    while (!st.empty()) {
//        TreeNode* node = st.top();
//        st.pop();
//        if (node != NULL) result.push_back(node->val);
//        else continue;
//        st.push(node->left); // 相对于前序遍历，这更改一下入栈顺序
//        st.push(node->right);
//    }
//    reverse(result.begin(), result.end()); // 将结果反转之后就是左右中的顺序了
//    return result;
//}

//int result;
//int n = 0;
//int kthSmallest(TreeNode* root, int k) {
//    if (root)
//    {
//        kthSmallest(root->left, k);
//        ++n;
//        if (k == n)
//        {
//            result = root->val;
//        }
//        kthSmallest(root->right, k);
//    }
//
//    return result;
//}

//int count(TreeNode* root)
////获取结点个数
//{
//    if (root == nullptr) return 0;
//    return 1 + count(root->left) + count(root->right);
//}
//int kthSmallest(TreeNode* root, int k) {
//    int cnt = count(root->left);//左子树结点个数
//    if (cnt == k - 1) return root->val;
//    else if (cnt > k - 1) return kthSmallest(root->left, k);
//    else return kthSmallest(root->right, k - cnt - 1);
//
//}

//TreeNode* invertTree(TreeNode* root) {
//    if (root != nullptr)
//    {
//        TreeNode* temp = root->left;
//        root->left = root->right;
//        root->right = temp;
//
//        invertTree(root->left);
//        invertTree(root->right);
//    }
//
//    return root;
//}

//vector<int> vis;
//void backtrack(vector<int>& nums, vector<vector<int>>& ans, int idx, vector<int>& perm) {
//    if (idx == nums.size()) {
//        ans.emplace_back(perm);
//        return;
//    }
//    for (int i = 0; i < (int)nums.size(); ++i) {
//        if (vis[i] || (i > 0 && nums[i] == nums[i - 1] && !vis[i - 1])) {
//            continue;
//        }
//        perm.emplace_back(nums[i]);
//        vis[i] = 1;
//        backtrack(nums, ans, idx + 1, perm);
//        vis[i] = 0;
//        perm.pop_back();
//    }
//}
//
//vector<vector<int>> permuteUnique(vector<int>& nums) {
//    vector<vector<int>> ans;
//    vector<int> perm;
//    vis.resize(nums.size());
//    sort(nums.begin(), nums.end());
//    backtrack(nums, ans, 0, perm);
//    return ans;
//}

//int sumOfLeftLeaves(TreeNode* root) {
//    int result = 0;
//    
//    if (root == nullptr)
//    {
//        return result;
//    }
//
//    queue<TreeNode*> que;
//    que.push(root);
//    while (!que.empty())
//    {
//        root = que.front();
//        que.pop();
//
//        if (root->left != nullptr)
//        {   
//            if (root->left->left == nullptr && root->left->right == nullptr)
//            {
//                result += root->left->val;
//            }
//            else
//            {
//                que.push(root->left);
//            }
//        }
//        if (root->right != nullptr)
//        {
//            que.push(root->right);
//        }
//    }
//
//    return result;
//}

//string intToRoman(int num) {
//	int values[] = { 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1 };
//	string reps[] = { "M", "CM", "D", "CD", "C", "XC", "L", "XL", "X", "IX", "V", "IV", "I" };
//
//	string res;
//	for (int i = 0; i < 13; i++)  //这里不使用图里的count了，一遍一遍来就行了
//		while (num >= values[i])
//		{
//			num -= values[i];
//			res += reps[i];
//		}
//	return res;
//}

vector<int> t;
vector<vector<int>> ans;

void dfs(int cur, vector<int>& nums) {
    if (cur == nums.size()) {
        ans.push_back(t);
        return;
    }
    t.push_back(nums[cur]);
    dfs(cur + 1, nums);
    t.pop_back();
    dfs(cur + 1, nums);
}

vector<vector<int>> subsets(vector<int>& nums) {
    dfs(0, nums);
    return ans;
}

int main()
{
    std::cout << "Hello World!\n";

    /*ector<int> a{ 0,1,2,2,3,0,4,2 };
    removeElement(a, 2);*/

    /*vector<int> a{ 1,1,2 };
    removeDuplicates(a);*/

    /*vector<int> a{ 2,3,6,7 };
    combinationSum(a, 8);*/

    /*combine(4, 3);*/

    /*cout << rangeBitwiseAnd(2147483646,2147483647) << endl;*/

    /*cout << countSubstrings("aaa") << endl;*/

    /*vector<int> a{ -10,-3,0,5,9 };
    sortedArrayToBST(a);*/

    /*vector<vector<int>> a{ {1,1,1} ,{1,1,0},{1,0,1} };
    floodFill(a, 1, 1, 2);*/

    //cout << multiply("12", "65") << endl;

    /*Node a1(1);
    Node a2(2);
    Node a3(3);
    Node a4(4);
    a1.neighbors.push_back(&a2);
    a1.neighbors.push_back(&a4);

    a2.neighbors.push_back(&a1);
    a2.neighbors.push_back(&a3);

    a3.neighbors.push_back(&a2);
    a3.neighbors.push_back(&a4);

    a4.neighbors.push_back(&a1);
    a4.neighbors.push_back(&a3);

    cloneGraph(&a1);*/

    //cout << countBinarySubstrings("00110011") << endl;

    /*restoreIpAddresses("25525511135");*/

    /*TreeNode a1 = 1;
    TreeNode a2 = 1;
    TreeNode a3 = 10;
    TreeNode a4 = 1;
    TreeNode a5 = 9;

    a1.left = nullptr;
    a1.right = &a2;

    a2.left = &a3;
    a2.right = &a4;

    a3.left = &a5;

    ListNode b1(1);
    ListNode b2(10);
    b1.next = &b2;

    cout << isSubPath(&b1, &a1) << endl;*/

    /*TreeNode a1 = 1;
    TreeNode a2 = 3;
    TreeNode a3 = 2;

    a1.left = &a2;
    a1.right = nullptr;

    a2.left = nullptr;
    a2.right = &a3;

    recoverTree(&a1);*/

    /*TreeNode a1 = 4;
    TreeNode a2 = 1;
    TreeNode a3 = 2;
    TreeNode a4 = 3;

    a1.left = &a2;
    a1.right = nullptr;

    a2.left = &a3;
    a2.right = nullptr;

    a3.left = &a4;

    cout << rob(&a1) << endl;*/

    /*vector<vector<int>> a{ {1,0 },{2,0 },{3,1 },{3,2} };
    findOrder(4, a);*/

    /*vector<vector<int>>a{ {1,0} };
    cout << canFinish(2,a) << endl;*/

    /*cout << addStrings("1","9") << endl;*/

    /*vector<int> a{ 0,0,0,0,0,0 };
    threeSum(a);*/

    //vector<vector<int>> a{ {0,1,0,0} ,{1,1,1,0},{0,1,0,0},{1,1,0,0} };
    //vector<vector<int>> a{ {0,1} };
    /*vector<vector<int>> a{ {0},{1}};
    cout << islandPerimeter(a) << endl;*/

    /*vector<int> a{5,6,5};
    pairSums(a, 11);*/

    /*vector<int>a{ 4,9,5 };
    vector<int>b{ 9,4,9,8,4 };
    intersect(a, b);*/

    /*vector<int>a{ 1,1,2,2,3,3 };
    cout << distributeCandies(a) << endl;*/

    /*solveNQueens(5);*/

    /*generateParenthesis(3);*/

    //vector<int> a{ 1,3,2,2,1 };
    /*vector<int> a{ 1,2,87,87,87,2,1 };
    cout << candy(a) << endl;*/

    /*vector<int> a{ 1,2,3 };
    vector<int> b{ 3 };
    cout << findContentChildren(a,b) << endl;*/

    /*vector<int> a{ 5, 1, 2, 3, 4 };
    vector<int> b{ 4, 4, 1, 5, 1 };
    cout << canCompleteCircuit(a,b) << endl;*/

    //cout << isSubsequence("axc", "aahbbgcdc") << endl;

	/*vector<string> d{ "looked", "just", "like", "her", "brother" };

	cout << respace(d, "jesslookedjustliketimherbrother") << endl;*/


    /*divingBoard(1,2,3);*/

    //vector<int> a{ 2,3,1,1,4 };
    //vector<int> a{ 2,0,2,0,1 };
    //vector<int> a{ 3, 4, 3, 2, 5, 4, 3 };
    //vector<int> a{ 1, 2, 1, 1, 1};
    /*cout << jump(a) << endl;*/

    //vector<int> a{ 3,2,1,0,4 };
    //vector<int> a{ 2,2,1,0,4 };
    //vector<int> a{ 2,3,1,1,4 };
    //vector<int> a{ 2,0};
    //vector<int> a{ 1,0,1,0};
    //vector<int> a{3,0,8,2,0,0,1};
    //cout << canJump(a) << endl;

    //vector<int> a{ 7,1,5,3,6,4 };
    //vector<int> a{ 7,6,4,3,1 };
    //cout << maxProfit(a) << endl;

    /*vector<int>a{ 5334,6299,4199,9663,8945,3566,9509,3124,6026,6250,7475,5420,9201,9501,38,5897,4411,6638,9845,161,9563,8854,3731,5564,5331,4294,3275,1972,1521,2377,3701,6462,6778,187,9778,758,550,7510,6225,8691,3666,4622,9722,8011,7247,575,5431,4777,4032,8682,5888,8047,3562,9462,6501,7855,505,4675,6973,493,1374,3227,1244,7364,2298,3244,8627,5102,6375,8653,1820,3857,7195,7830,4461,7821,5037,2918,4279,2791,1500,9858,6915,5156,970,1471,5296,1688,578,7266,4182,1430,4985,5730,7941,3880,607,8776,1348,2974,1094,6733,5177,4975,5421,8190,8255,9112,8651,2797,335,8677,3754,893,1818,8479,5875,1695,8295,7993,7037,8546,7906,4102,7279,1407,2462,4425,2148,2925,3903,5447,5893,3534,3663,8307,8679,8474,1202,3474,2961,1149,7451,4279,7875,5692,6186,8109,7763,7798,2250,2969,7974,9781,7741,4914,5446,1861,8914,2544,5683,8952,6745,4870,1848,7887,6448,7873,128,3281,794,1965,7036,8094,1211,9450,6981,4244,2418,8610,8681,2402,2904,7712,3252,5029,3004,5526,6965,8866,2764,600,631,9075,2631,3411,2737,2328,652,494,6556,9391,4517,8934,8892,4561,9331,1386,4636,9627,5435,9272,110,413,9706,5470,5008,1706,7045,9648,7505,6968,7509,3120,7869,6776,6434,7994,5441,288,492,1617,3274,7019,5575,6664,6056,7069,1996,9581,3103,9266,2554,7471,4251,4320,4749,649,2617,3018,4332,415,2243,1924,69,5902,3602,2925,6542,345,4657,9034,8977,6799,8397,1187,3678,4921,6518,851,6941,6920,259,4503,2637,7438,3893,5042,8552,6661,5043,9555,9095,4123,142,1446,8047,6234,1199,8848,5656,1910,3430,2843,8043,9156,7838,2332,9634,2410,2958,3431,4270,1420,4227,7712,6648,1607,1575,3741,1493,7770,3018,5398,6215,8601,6244,7551,2587,2254,3607,1147,5184,9173,8680,8610,1597,1763,7914,3441,7006,1318,7044,7267,8206,9684,4814,9748,4497,2239};
    cout << splitArray(a,9) << endl;*/

    //cout << divisorGame(3) << endl;

    //vector<vector<int>>a{ {1,3,1},{1,5,1},{4,2,1} };
    //vector<vector<int>>a{ {1,2},{1,1} };
    /*vector<vector<int>>a{ {1,2,5},{3,2,1}};
    cout << minPathSum(a) << endl;*/

    //vector<int> a = { 10,1,10,10,10 };
    /*vector<int> a = { 10,10,10,1,10};
    cout << minArray(a) << endl;*/

    /*generateTrees(5);*/

    /*vector<int> a = { 2,3,4};
    twoSum(a, 6);*/


    /*vector<int> a = { 3,1,5,8 };
    cout << maxCoins(a) << endl;*/

    /*TreeNode a1 = 1;
    TreeNode a2 = 2;
    TreeNode a3 = 3;

    a1.left = nullptr;
    a1.right = &a2;

    a2.left = &a3;
    a2.right = nullptr;


    inorderTraversal(&a1);*/

    /*TreeNode a1 = 5;
    TreeNode a2 = 1;
    TreeNode a3 = 4;
    TreeNode a4 = 3;
    TreeNode a5 = 6;

    a1.left = &a2;
    a1.right = &a3;

    a2.left = nullptr;
    a2.right = nullptr;

    a3.left = &a4;
    a3.right = &a5;

    cout << isValidBST(&a1) << endl;*/

    /*TreeNode a1 = 1;
    TreeNode a2 = 2;
    a1.left = &a2;
    TreeNode a3 = 1;
    a2.right = &a3;

    TreeNode b1 = 1;
    TreeNode b2 = 1;
    b1.left = &b2;
    TreeNode b3 = 2;
    b2.right = &b3;*/
    //cout << isSameTree(&a1, &b1) << endl;


    //vector<int> a = { 1,3 };
    //vector<int> a = {1,3,5,6};
    //cout << searchInsert(a,2) << endl;

    /*cout << numTrees(3);*/

    //vector<vector<int>> a = { {2},{3,4},{6,5,7},{4,1,8,3} };
    //vector<vector<int>> a = { {-1},{-2,-3}};
    /*vector<vector<int>> a = { {-10} };
    cout << minimumTotal(a) << endl;*/

    /*cout << reverse(1563847412);*/

    /*cout << convert("AB",1) << endl;*/

    /*vector<int> a = { 1, 5, 11, 5 };
    bool b = canPartition(a);
        cout << b << endl;*/

    /*int a = integerBreak(10);
    cout << a << endl;*/

    /*int a = numSquares(12);*/

    //cout << longestPalindrome("cbbd") << endl;

    //int length = lengthOfLongestSubstring("abcabcbb");
    //int length = lengthOfLongestSubstring("bbbbb");    
    //int length = lengthOfLongestSubstring("pwwkew");
    //int length = lengthOfLongestSubstring("cdd");
    //int length = lengthOfLongestSubstring("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789!\"#$%&'()*+,-./:;<=>?@[\\]^_`{|}~");
    //cout << length << endl;


    /*ListNode a1 = 2;
    
    ListNode a2 = 4;
    a1.next = &a2;
    ListNode a3 = 3;
    a2.next = &a3;
    ListNode b1 = 5;
    ListNode b2 = 6;
    b1.next = &b2;
    ListNode b3 = 4;
    b2.next = &b3;

    ListNode* result = addTwoNumbers(&a1, &b1);*/

    cout << endl;
}



// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
