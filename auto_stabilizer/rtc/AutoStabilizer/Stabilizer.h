#ifndef Stabilizer_H
#define Stabilizer_H

#include "GaitParam.h"
#include <prioritized_qp_osqp/prioritized_qp_osqp.h>
#include <cnoid/JointPath>

class Stabilizer{
public:
  // Stabilizerでしか使わないパラメータ
  std::vector<double> bodyAttitudeControlGain=std::vector<double>{0.5, 0.5}; // 要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[/s]. 0以上
  std::vector<double> bodyAttitudeControlTimeConst=std::vector<double>{1000, 1000}; // 要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[s]. 0より大きい
  std::vector<double> bodyAttitudeControlCompensationLimit=std::vector<double>{0.7,0.7}; //要素数2. gaitParam.footMidCoords座標系X軸, Y軸. 単位は[rad]. STが動いている間は変更されない. 0以上

  std::vector<std::vector<double> > supportPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > supportDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > landingPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > landingDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > swingPgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  std::vector<std::vector<double> > swingDgain = std::vector<std::vector<double> >(2); // 要素数2. [rleg, lleg]. rootLinkから各endeffectorまでの各関節のゲイン. 0~100
  double swing2LandingTransitionTime = 0.05; // [s]. 0より大きい
  double landing2SupportTransitionTime = 0.1; // [s]. 0より大きい
  double support2SwingTransitionTime = 0.2; // [s]. 0より大きい

  std::vector<cnoid::Vector6> ee_K; // 要素数EndEffectors. EndEffector frame. endEffector origin. 0以上
  std::vector<cnoid::Vector6> ee_D; // 要素数EndEffectors. EndEffector frame. endEffector origin. 0以上
  cnoid::Vector3 com_K; // generate frame. 
  cnoid::Vector3 com_D; // generate frame.
  cnoid::VectorXd refAngle_K; // generate frame.
  cnoid::VectorXd refAngle_D; // generate frame.
  double ee_dv_limit = 20.0; // 分解加速度制御でのタスク空間でのフィードバック込みの加速度ノルム上限
  double ee_dw_limit = 20.0; // 分解加速度制御でのタスク空間でのフィードバック込みの角加速度ノルム上限
  double defaultDdqLimit = 50; // qpにいれる全関節共通のデフォルト関節角加速度リミット. 特にrootについてはこれを使う
  std::vector<double> ddq_limit; // qpに入れる関節角加速度リミット. 電流リミット*トルク定数*ギア比= トルクリミットを関節イナーシャで割った値. 要素数はnumJointsなのでrootを含まない. 大きすぎて(低くて400m/s^2)実際は機能していない
  std::vector<double> torque_limit; // u に入れる関節トルクリミット. 電流リミット*トルク定数*ギア比
  std::vector<double> torque_ratcheting_limit = {220, 450, 1000, 1000, 220, 220, // 右足
                                                 220, 450, 1000, 1000, 220, 220, // 左足
                                                 1000, 1000, 220, 220, 220, // 胴体
                                                 220, 220, 220, 220, 220, 220, 220, 220, // 右腕
                                                 220, 220, 220, 220, 220, 220, 220, 220, // 左腕
                                                 10, 10, 10, 10, 10, 10, // 右指
                                                 10, 10, 10, 10, 10, 10, // 左指
                                                 };
  void init(const GaitParam& gaitParam, cnoid::BodyPtr& actRobotTqc){
    for(int i=0;i<NUM_LEGS;i++){
      cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
      if(jointPath.numJoints() == 6){
        /*supportPgain[i] = {5,15,10,5,0.2,0.2};
        supportDgain[i] = {10,20,20,10,5,5};
        landingPgain[i] = {5,15,1,1,0.2,0.2};
        landingDgain[i] = {10,10,10,10,5,5};
        swingPgain[i] = {5,30,20,10,5,5};
        swingDgain[i] = {10,30,20,20,30,30};*/
        supportPgain[i] = {0,0,0,0,0,0};
        supportDgain[i] = {0,0,0,0,0,0};
        landingPgain[i] = {0,0,0,0,0,0};
        landingDgain[i] = {0,0,0,0,0,0};
        swingPgain[i] = {0,0,0,0,0,0};
        swingDgain[i] = {0,0,0,0,0,0};
	// Yamamoto param
	// JAXON_RED固有（モータードライバ種類依存）
	// pgain = {100,100,100,300,100,100}
	// dgain = {20,2,2,5,3,3}
	/*supportPgain[i] = {0.31,0.12,0.15,0.46,0.21,0.30};
        supportDgain[i] = {1.693,0.066,0.085,0.212,0.176,0.251};
        landingPgain[i] = {0.31,0.12,0.15,0.46,0.21,0.30};
        landingDgain[i] = {1.693,0.066,0.085,0.212,0.176,0.251};
        swingPgain[i] = {0.31,0.12,0.15,0.46,0.21,0.30};
        swingDgain[i] = {1.693,0.066,0.085,0.212,0.176,0.251};*/
      }else{
        supportPgain[i].resize(jointPath.numJoints(), 100.0);
        supportDgain[i].resize(jointPath.numJoints(), 100.0);
        landingPgain[i].resize(jointPath.numJoints(), 100.0);
        landingDgain[i].resize(jointPath.numJoints(), 100.0);
        swingPgain[i].resize(jointPath.numJoints(), 100.0);
        swingDgain[i].resize(jointPath.numJoints(), 100.0);
      }
    }
    this->ddq_limit.resize(gaitParam.actRobotTqc->numJoints());
    this->torque_limit.resize(gaitParam.actRobotTqc->numJoints());
    for (int i=0;i<gaitParam.actRobotTqc->numJoints();i++){
      double climit = gaitParam.actRobotTqc->joint(i)->info<double>("climit");
      if(climit >= 200) climit = 10; // 大きすぎるor設定されていないので適当
      double gearRatio = gaitParam.actRobotTqc->joint(i)->info<double>("gearRatio");
      double torqueConst = gaitParam.actRobotTqc->joint(i)->info<double>("torqueConst");
      cnoid::Matrix3 Inertia = gaitParam.actRobotTqc->joint(i)->I();
      cnoid::Vector3 axis = gaitParam.actRobotTqc->joint(i)->jointAxis();
      this->torque_limit[i] = std::max(climit * gearRatio * torqueConst, this->torque_ratcheting_limit[i]);
      this->ddq_limit[i] = climit * gearRatio * torqueConst / (Inertia * axis).norm();
    }
    for(int i=0;i<gaitParam.eeName.size();i++){
      cnoid::Vector6 defaultEED;
      if(i<NUM_LEGS){
	defaultEED << 30, 30, 30, 20, 20, 20;
      }else{
	defaultEED << 10, 10, 10, 10, 10, 10;
      }
      this->ee_D.push_back(defaultEED);
      cnoid::Vector6 defaultEEK;
      if(i<NUM_LEGS){
	defaultEEK << 200, 200, 200, 100, 100, 100;
      }else{
	defaultEEK << 50, 50, 50, 20, 20, 20;
      }
      this->ee_K.push_back(defaultEEK);
    }
    cnoid::Vector3 defaultComK;
    defaultComK << 0, 0, 0;
    this->com_K = defaultComK;
    cnoid::Vector3 defaultComD;
    defaultComD << 0, 0, 0;
    this->com_D = defaultComD;
    this->refAngle_K = cnoid::VectorXd::Zero(6 + gaitParam.actRobotTqc->numJoints());
    this->refAngle_D = cnoid::VectorXd::Zero(6 + gaitParam.actRobotTqc->numJoints());
    cnoid::Vector6 defaultRootK;
    defaultRootK << 10, 10, 10, 100, 100, 100;
    this->refAngle_K.head<6>() = defaultRootK;
    for (int i=0;i<gaitParam.actRobotTqc->numJoints();i++){
      this->refAngle_K[6+i] = 1;
      if((i==12) || (i==13) || (i==14)) this->refAngle_K[6+i] = 100; // 腰roll pitch yaw
    }
    cnoid::Vector6 defaultRootD;
    defaultRootD << 1, 1, 1, 10, 10, 10;
    this->refAngle_D.head<6>() = defaultRootD;
    for (int i=0;i<gaitParam.actRobotTqc->numJoints();i++){
      this->refAngle_D[6+i] = 1;
      if((i==12) || (i==13) || (i==14)) this->refAngle_D[6+i] = 10; // 腰roll pitch yaw
    }
  }
protected:
  // 計算高速化のためのキャッシュ. クリアしなくても別に副作用はない.
  // for calcWrench
  mutable std::shared_ptr<prioritized_qp_osqp::Task> constraintTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> tgtZmpTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> copTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  // for calcTorque
  mutable std::shared_ptr<prioritized_qp_osqp::Task> eeComTask_ = std::make_shared<prioritized_qp_osqp::Task>();
  mutable std::shared_ptr<prioritized_qp_osqp::Task> jointPDTask_ = std::make_shared<prioritized_qp_osqp::Task>();
public:
  void initStabilizerOutput(const GaitParam& gaitParam,
                            cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Vector3& o_stTargetZmp, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage) const;

  bool execStabilizer(const GaitParam& gaitParam, double dt, bool useActState,
                      cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const;

  bool calcResolvedAccelationControl(const GaitParam& gaitParam, double dt, bool useActState, cnoid::BodyPtr& actRobotTqc, 
				     cnoid::Vector3& o_stTargetZmp, std::vector<cnoid::Vector6>& o_stEETargetWrench,
				     std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage,
				     Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd, cnoid::Vector6& prev_rootd) const;

protected:
  bool moveBasePosRotForBodyRPYControl(double dt, const GaitParam& gaitParam, bool useActState,
                                       cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const;
  bool calcZMP(const GaitParam& gaitParam, double dt, bool useActState,
               cnoid::Vector3& o_tgtZmp/*generate座標系*/, cnoid::Vector3& o_tgtForce/*generate座標系*/, cnoid::Vector3& o_tgtCogAcc /*generate座標系*/) const;
  bool calcWrench(const GaitParam& gaitParam, const cnoid::Vector3& tgtZmp/*generate座標系*/, const cnoid::Vector3& tgtForce/*generate座標系. ロボットが受ける力*/, bool useActState, cnoid::BodyPtr& actRobotTqc, 
                  std::vector<cnoid::Vector6>& o_tgtEEWrench /* 要素数EndEffector数. generate座標系. EndEffector origin*/) const;
  bool calcTorque(double dt, const GaitParam& gaitParam, bool useActState, cnoid::BodyPtr& actRobotTqc, const cnoid::Vector3& targetCogAcc, 
                  std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage,
		  cnoid::Vector3& root2CogForce, Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd, cnoid::Vector6& prev_rootd) const;

};

#endif
