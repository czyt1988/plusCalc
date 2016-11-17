function st = getStraightPipePulsationCalcDefaultStrct()
%获取直管气流脉动计算的默认输入结构体
%   
st.massFlowE = nan;
st.k = nan;
st.a = nan;
st.S = nan;
st.Dpipe = nan;
st.isDamping = 0;
st.coeffFriction = nan;
st.meanFlowVelocity = nan;
st.mach = nan;
st.notMach = 0;
st.frequency = nan;
st.time = nan;
st.L = nan;
st.sectionL = nan;
end

