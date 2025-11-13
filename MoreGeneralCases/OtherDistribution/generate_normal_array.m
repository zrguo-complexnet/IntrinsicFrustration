function result = generate_normal_array(N,center,sigma)
    % 正态分布数组生成（标准正态分布，均值0，标准差1）
    result = norminv((1:(N+1)) / (N+1), center, sigma);
    result(N+1)=[];
end
