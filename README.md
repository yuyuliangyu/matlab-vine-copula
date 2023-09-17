# data
这是一个用matlab编写的vine-copula
参考的文献是
Selecting and estimating regular vine copulae and application to financial returns(J. Dißmann, E. C. Brechmann∗, C. Czado, D. Kurowicka)
提供的数据是亚洲某五个地区的指数，仅供参考
#
需要注意的是：
成功设计了vine-copula的拟合函数，但是并没有编写vine-copula的模拟函数。
基本的分布函数只涉及正太，t，Clayton等基础分布
有相关领域的人员编写了用matlab结合c++的vine-copula代码，而此文件全是matlab代码。
#
其实以上注意事项原因归结于matlab并没有成熟的统计分布相关的函数（但是c++有相应的强大工具），这对于求解vine-copula中的逆分布是一个极大的困难。而此文件只设计了一些基础的分布函数，例如正太，t，clayton等，
这些都是在matlab2018中有的数据。
#
希望相关学者专家可以提供思路，欢迎指证~
by yuyu , a student from zhejiang gongshang university
